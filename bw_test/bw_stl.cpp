
#include <iostream>
#include <time.h>
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




double Prpoiss(int cn,  int cov, int Hmean) {
    double result=0;
    const double epsi=1e-99;
    if(cov > std::floor(30.497*Hmean) ){
        if(cn!= 31)
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
            if (emisP[i][j]==0)
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
void train(vector<size_t> &observations,
           vector<double> &startP,
           vector<vector<double>> &transP,
           vector<vector<double>> &emisP,
           const int max_iter )
{

    // 1. Initialization
    int iters = 0;
   // int T = observations.cols; // number of element per sequence
   // int C = observations.rows; // number of sequences
    int nStates = transP.rows; // number of states | also nStates = transP.cols | transP = A = {aij} - NxN
    int nObservations = emisP.cols; // number of observations | emisP = B = {bj(k)} - NxM
    
    correctModel(transP,emisP,startP);


    //cv::Mat FTRANS,FINIT,FEMIS;
    vector<vector<double>> FTRANS(nStates, vector<double>(nStates));
    FTRANS = transP;

    vector<vector<double>> FEMIS(nStates, vector<double>(max_obs  ));
    FEMIS = emisP;

    vector<double> FINIT(nStates);
    FINIT = startP;
    

    //-1 * (std::numeric_limits<double>::max());
    double logProb = -DBL_MAX;
    double oldLogProb;
    int data = 0;
    do {
        oldLogProb = logProb;

        // compute a0
        vector<vector<double>> a(nStates, vector<double>(nObservations));
        vector<double> c(nObservations);
        c[0]=0;
        
//        cv::Mat a(nStates,T,CV_64F);
  //      cv::Mat c(1,T,CV_64F); c.at<double>(0,0) = 0;
   
        for (int i=0;i<nStates;i++)
        {
            a[i][0] = startP[0][i]*emisP[i][observations[0]];
            c[0] += a[i][0];
        }
        // scale the a0(i)
        
        c.at<double>(0,0) = 1/c.at<double>(0,0);
        for (int i=0;i<nStates;i++)
            a[i][0] *= c[0];

        // 2. The a-pass
        // compute at(i)
        for (int t=1;t<T;t++)
        {
            c[t] = 0;
            for (int i=0;i<nStates;i++)
            {
                a[i][t] = 0;
                for (int j=0;j<nStates;j++)
                    a[i][t] += a[i][t-1]*transP[j][i];
                a[i][t] = a[i][t] * emisP[i][observations[t]];
                c[t] += a[i][t];
            }
            // scale at(i)
            c[t] = 1/c[t];
            for (int i=0;i<nStates;i++)
                a[i][t]=c[t] * a[i][t];
        }
        // 3. The B-pass
        vector<vector<double>> b(nStates, vector<double>(nObservations));
   //     cv::Mat b(nStates,T,CV_64F);
        // Let Bt-1(i) = 1 scaled by Ct-1
        for (int i=0;i<nStates;i++)
            b[i][T-1] = c[T-1];
        // B-pass
        for (int t=T-2;t>-1;t--)
            for (int i=0;i<nStates;i++)
            {
                b[i][t] = 0;
                for (int j=0;j<nStates;j++)
                    b[i][t] += transP[i][j] * emisP[j][observations[t+1]] * b[j][t+1];
                // scale Bt(i) with same scale factor as at(i)
                b[i][t] *= c[t];
            }
        // 4. Compute  Yt(i,j) and Yt(i)
        double denom;
        int index;
        vector<vector<double>> YN(nStates, vector<double>(nObservations));
        vector<vector<double>> YNN(nStates*nStates, vector<double>(nObservations));
//        cv::Mat YN(nStates,T,CV_64F);
  //      cv::Mat YNN(nStates*nStates,T,CV_64F);
        for (int t=0;t<T-1;t++)
        {
            denom = 0;
            for (int i=0;i<nStates;i++)
                for (int j=0;j<nStates;j++)
                    denom += a[i][t] * transP[i][j] * emisP[j][observations[t+1]] * b[j][t+1];
            index = 0;
            for (int i=0;i<nStates;i++)
            {
                YN[i][t] = 0;
                for (int j=0;j<nStates;j++)
                {
                    YNN[index][t] = (a[i][t] * transP[i][j] * emisP[j][observations[t+1]] * b[j][t+1] ) / denom;
                    YN[i][t]+=YNN[index][t];
                    index++;
                }
            }
        }
        // 5. Re-estimate A,B and pi
        // re-estimate pi
        for (int i=0;i<nStates;i++)
            startP[0][i] = YN[i][0];
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
                    numer += YNN[index][t];
                    denom += YN[i][t];
                }
                transP[i][j] = numer/denom;
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
                    if (observations[t] == j)
                        numer+= YN[i][t];
                    denom += YN[i][t];
                }
                emisP[i][j] = numer/denom;
            }
        correctModel(transP,emisP,startP);
        
        /*
        FTRANS = (FTRANS*(data+1)+transP)/(data+2);
        FEMIS = (FEMIS*(data+1)+emisP)/(data+2);
        FINIT = (FINIT*(data+1)+startP)/(data+2);
        */
        
        
        // 6. Compute log[P(O|y)]
        logProb = 0;
        for (int i=0;i<T;i++)
            logProb += log(c[i]);
        logProb *= -1;
        
        /* 7. To iterate or not
        data++;
        if (data >= C)
        {
            data = 0;
            iters++;
        }
        */
        
        //FIX ITER AND LOGPROB DIFFERENCE
        
    } while (iters<max_iter && logProb>oldLogProb);
    correctModel(FTRANS,FEMIS,FINIT);
    
    vector<vector<double>> FtransP(nStates, vector<double>(nStates));
    FtransP = FTRANS;

    vector<vector<double>> FemisP(nStates, vector<double>(max_obs  ));
    FemisP = FEMIS;

    vector<double> FstartP(nStates);
    FstartP = FINIT;
}//train

//FIX :NEED TO RETURN FxxxP



int main(int argc, const char * argv[])
{
    
    if (argc != 8) {
        std::cerr << "usage: " << argv[0] << " <observation> <nStates> <Haploid Mean> <output prefix> <max_cov> <max_iter> <delta>" << endl;
        return EXIT_FAILURE;
    }
    
    const double epsi=1e-99;

    const string prefix_file(argv[4]);
    
    
    const size_t mean=std::stoi(argv[3]);

    //number of states
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
    const string filename2(argv[1]);
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
    size_t nObservations=observations.size();

    
    
    //max cov value observed or upper cov bound
    size_t max_obs = std::min(std::ceil(30.497*mean),std::stoi(argv[5]));
    
    //transition probabilities
    vector<vector<double>> transP(nStates, vector<double>(nStates));
    for (size_t i=0;i<nStates;i++)
    {
        for (size_t j=0;j<nStates;j++)
        {
            if(i==j)
            {
                transP[i][j]=1 - ( ( ((double) (nStates)) -1) * (epsi) );
            //        cout<<" "<<transP[i][j];
            }
        else
        {
            transP[i][j]=epsi;
            //        cout<<" "<<transP[i][j];
        }
        }//cout<<endl;
    }
    
    //emission probabilities
    vector<vector<double>> emisP(nStates, vector<double>(max_obs  ));
    for (size_t i=0;i<nStates;i++)
    {
        for (size_t j=0;j<=max_obs;j++)
        {
            emisP[i][j]=Prpoiss( (int) i ,j, (int) mean );
        }//cout<<endl;
    }




    vector<double> startP(nStates);
    for(size_t i=0;i<(nStates);i++)
    {
        startP[i]=1./(nStates);
    }
    
    const int max_it=std::stoi(argv[6]);
    train(observations,transP, emisP , startP, max_it);
    printModel(TRGUESS,EMITGUESS,INITGUESS);

    std::cout << "\ndone.\n";

    return 0;
}
