

#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <cassert>
#include <cmath>
#include <random>
#include <boost/math/distributions/poisson.hpp>
using boost::math::poisson;
using boost::math::pdf;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;

#include <iostream>
#include <time.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>



double Prpoiss(int cn,  int cov, int Hmean) {
    double result=0;
    const double epsi=1e-99;
    if(cov > 20.49*Hmean){
        if(cn!= 21)
            result=epsi;
        else
            result=1-epsi;
    }
    else if(cn!=0){
        poisson distribution(cn*Hmean);
        result=pdf(distribution, cov);
    }
    else{
        poisson distribution(0.1*Hmean);
        result=pdf(distribution, cov);
    }
    return result;
}



static double emissionPr( size_t index, size_t state, vector<size_t> & observations,size_t mean  ){

    int cov = observations[index];

    /*
     int cov = coverage[index];

     if (as.integer( obs[2])  ==1)
     {
     if ( stt[2]  =="1")
     result=Prpoiss(  as.integer(stt[1]) ,cov,30 ) * Prclip(3,3);
     else
     result=Prpoiss( as.integer(stt[1]) ,cov,30  ) * Prclip(0,3);
     }
     else
     {
     if ( stt[2]  =="1")
     result=Prpoiss(  as.integer(stt[1]) ,coverage,30  ) * Prclip(0,3);
     else
     result=Prpoiss( as.integer(stt[1]) ,coverage,30  ) * Prclip(3,3) ;
     }
     */

    double result=Prpoiss( (int) state ,cov, (int) mean );
    return result;
}

static void correctModel(cv::Mat &transP, cv::Mat &emisP, cv::Mat &startP)
{
    double eps = 1e-30;
    for (int i=0;i<emisP.rows;i++)
        for (int j=0;j<emisP.cols;j++)
            if (emisP.at<double>(i,j)==0)
                emisP.at<double>(i,j)=eps;
    
    for (int i=0;i<transP.rows;i++)
        for (int j=0;j<transP.cols;j++)
            if (transP.at<double>(i,j)==0)
                transP.at<double>(i,j)=eps;
    
    for (int i=0;i<startP.cols;i++)
        if (startP.at<double>(0,i)==0)
            startP.at<double>(0,i)=eps;
    
    double sum;
    for (int i=0;i<transP.rows;i++)
    {
        sum = 0;
        for (int j=0;j<transP.cols;j++)
            sum+=transP.at<double>(i,j);
        for (int j=0;j<transP.cols;j++)
            transP.at<double>(i,j)/=sum;
    }
    for (int i=0;i<emisP.rows;i++)
    {
        sum = 0;
        for (int j=0;j<emisP.cols;j++)
            sum+=emisP.at<double>(i,j);
        for (int j=0;j<emisP.cols;j++)
            emisP.at<double>(i,j)/=sum;
    }
    sum = 0;
    for (int j=0;j<startP.cols;j++)
        sum+=startP.at<double>(0,j);
    for (int j=0;j<startP.cols;j++)
        startP.at<double>(0,j)/=sum;
}//correctModel



static void printModel(const cv::Mat &transP,const cv::Mat &emisP,const cv::Mat &startP)
{
    std::cout << "\nTRANS: \n";
    for (int r=0;r<transP.rows;r++)
    {
        for (int c=0;c<transP.cols;c++)
            std::cout << transP.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\nEMIS: \n";
    for (int r=0;r<emisP.rows;r++)
    {
        for (int c=0;c<emisP.cols;c++)
            std::cout << emisP.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\nINIT: \n";
    for (int r=0;r<startP.rows;r++)
    {
        for (int c=0;c<startP.cols;c++)
            std::cout << startP.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}//printModel



/* Calculates maximum likelihood estimates of transition and emission probabilities from a sequence of emissions */
void train(    vector<size_t> &observations,vector<double> &startP,
           vector<vector<double>> &transP,   vector<vector<double>> &emisP, const int max_iter )
{

    // 1. Initialization
    int iters = 0;
   // int T = observations.cols; // number of element per sequence
   // int C = observations.rows; // number of sequences
    int nStates = transP.rows; // number of states | also nStates = transP.cols | transP = A = {aij} - NxN
    int nObservations = emisP.cols; // number of observations | emisP = B = {bj(k)} - NxM
    
    correctModel(transP,emisP,startP);


    cv::Mat FTRANS,FINIT,FEMIS;
        FTRANS = transP.clone();
        FEMIS = emisP.clone();
        FINIT = startP.clone();
    


    double logProb = -DBL_MAX;
    double oldLogProb;
    int data = 0;
    do {
        oldLogProb = logProb;

        // compute a0
        cv::Mat a(nStates,T,CV_64F);
        cv::Mat c(1,T,CV_64F); c.at<double>(0,0) = 0;
        for (int i=0;i<nStates;i++)
        {
            a.at<double>(i,0) = startP.at<double>(0,i)*emisP.at<double>(i,observations.at<int>(0,0));
            c.at<double>(0,0) += a.at<double>(i,0);
        }
        // scale the a0(i)
        c.at<double>(0,0) = 1/c.at<double>(0,0);
        for (int i=0;i<nStates;i++)
            a.at<double>(i,0) *= c.at<double>(0,0);

        // 2. The a-pass
        // compute at(i)
        for (int t=1;t<T;t++)
        {
            c.at<double>(0,t) = 0;
            for (int i=0;i<nStates;i++)
            {
                a.at<double>(i,t) = 0;
                for (int j=0;j<nStates;j++)
                    a.at<double>(i,t) += a.at<double>(i,t-1)*transP.at<double>(j,i);
                a.at<double>(i,t) = a.at<double>(i,t) * emisP.at<double>(i,observations.at<int>(data,t));
                c.at<double>(0,t)+=a.at<double>(i,t);
            }
            // scale at(i)
            c.at<double>(0,t) = 1/c.at<double>(0,t);
            for (int i=0;i<nStates;i++)
                a.at<double>(i,t)=c.at<double>(0,t)*a.at<double>(i,t);
        }
        // 3. The B-pass
        cv::Mat b(nStates,T,CV_64F);
        // Let Bt-1(i) = 1 scaled by Ct-1
        for (int i=0;i<nStates;i++)
            b.at<double>(i,T-1) = c.at<double>(0,T-1);
        // B-pass
        for (int t=T-2;t>-1;t--)
            for (int i=0;i<nStates;i++)
            {
                b.at<double>(i,t) = 0;
                for (int j=0;j<nStates;j++)
                    b.at<double>(i,t) += transP.at<double>(i,j)*emisP.at<double>(j,observations.at<int>(data,t+1))*b.at<double>(j,t+1);
                // scale Bt(i) with same scale factor as at(i)
                b.at<double>(i,t) *= c.at<double>(0,t);
            }
        // 4. Compute  Yt(i,j) and Yt(i)
        double denom;
        int index;
        cv::Mat YN(nStates,T,CV_64F);
        cv::Mat YNN(nStates*nStates,T,CV_64F);
        for (int t=0;t<T-1;t++)
        {
            denom = 0;
            for (int i=0;i<nStates;i++)
                for (int j=0;j<nStates;j++)
                    denom += a.at<double>(i,t)*transP.at<double>(i,j)*emisP.at<double>(j,observations.at<int>(data,t+1))*b.at<double>(j,t+1);
            index = 0;
            for (int i=0;i<nStates;i++)
            {
                YN.at<double>(i,t) = 0;
                for (int j=0;j<nStates;j++)
                {
                    YNN.at<double>(index,t) = (a.at<double>(i,t)*transP.at<double>(i,j)*emisP.at<double>(j,observations.at<int>(data,t+1))*b.at<double>(j,t+1))/denom;
                    YN.at<double>(i,t)+=YNN.at<double>(index,t);
                    index++;
                }
            }
        }
        // 5. Re-estimate A,B and pi
        // re-estimate pi
        for (int i=0;i<nStates;i++)
            startP.at<double>(0,i) = YN.at<double>(i,0);
        // re-estimate A
        double numer;
        index = 0;
        for (int i=0;i<nStates;i++)
            for (int j=0;j<nStates;j++)
            {
                numer = 0;
                denom = 0;
                for (int t=0;t<T-1;t++)
                {
                    numer += YNN.at<double>(index,t);
                    denom += YN.at<double>(i,t);
                }
                transP.at<double>(i,j) = numer/denom;
                index++;
            }
        // re-estimate B
        for (int i=0;i<nStates;i++)
            for (int j=0;j<nObservations;j++)
            {
                numer = 0;
                denom = 0;
                for (int t=0;t<T-1;t++)
                {
                    if (observations.at<int>(data,t)==j)
                        numer+=YN.at<double>(i,t);
                    denom += YN.at<double>(i,t);
                }
                emisP.at<double>(i,j) = numer/denom;
            }
        correctModel(transP,emisP,startP);
        FTRANS = (FTRANS*(data+1)+transP)/(data+2);
        FEMIS = (FEMIS*(data+1)+emisP)/(data+2);
        FINIT = (FINIT*(data+1)+startP)/(data+2);
        // 6. Compute log[P(O|y)]
        logProb = 0;
        for (int i=0;i<T;i++)
            logProb += log(c.at<double>(0,i));
        logProb *= -1;
        // 7. To iterate or not
        data++;
        if (data >= C)
        {
            data = 0;
            iters++;
        }
    } while (iters<max_iter && logProb>oldLogProb);
    correctModel(FTRANS,FEMIS,FINIT);
    transP = FTRANS.clone();
    emisP = FEMIS.clone();
    startP = FINIT.clone();
}//train





int main(int argc, const char * argv[])
{
    const size_t mean=std::stoi(argv[3]);
    
    
    
    if (argc != 5) {
        std::cerr << "usage: " << argv[0] << " <observation> <nStates> <Haploid Mean> <output prefix>" << endl;
        return EXIT_FAILURE;
    }
    
    
    
    const string filename2(argv[1]);
    const string prefix_file(argv[4]);
    const size_t nstates=ceil(std::stof(argv[2]));
    size_t nStates=0;
    if (nstates>30)
        nStates=31;
    else if (nStates<1)
        std::cerr << "nstates>0"<<endl;
    else
        nStates= nstates;
    
    //observations
    vector<size_t> observations;
size_t nObservations=observations.size();
    
    
    std::ifstream file;
    size_t inputString;
    file.open(filename2);
    int i=0;
    while(file >> inputString){
        observations.push_back(inputString);
        //observations[i]=inputString;
        i++;
    }
    file.close();

size_t max_obs=std::ceil(30.497*mean)
    
    
vector<vector<double>> transP(nStates, vector<double>(nStates));
for (size_t i=0;i<nStates;i++){
    for (size_t j=0;j<nStates;j++){
        if(i==j){
            transP[i][j]=1 - ( ( ((double) (nStates)) -1) * (epsi) );
            //        cout<<" "<<transP[i][j];
        }
        else{
            transP[i][j]=epsi;
            //        cout<<" "<<transP[i][j];
        }
    }//cout<<endl;
}
    
    
    vector<vector<double>> emisP(nStates, vector<double>(max_obs  ));
for (size_t i=0;i<nStates;i++){
  for (size_t j=0;j<distinct_nObservations;j++){
      emisP[i][j]=Prpoiss( (int) i ,j, (int) mean )
  }//cout<<endl;
}




const double epsi=1e-99;
vector<double> startP(nStates);

for(size_t i=0;i<(nStates);i++){
     startP[i]=1./(nStates);
}

    train(observations,100,transP, emisP , startP);
    printModel(TRGUESS,EMITGUESS,INITGUESS);




    std::cout << "\ndone.\n";




    return 0;
}
