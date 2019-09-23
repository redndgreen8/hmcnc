
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
        poisson distribution(0.01*Hmean);
        result=pdf(distribution, cov);
    }
    return result;
}



static double emissionPr( size_t index, size_t state, vector<size_t> & observations,size_t mean  ){

    int cov = observations[index];
    double result=Prpoiss( (int) state ,cov, (int) mean );
    return result;
}

static void correctModel(vector<vector<double>> &transP,
                         vector<vector<double>> &emisP,
                         vector<double> &startP,
                         int nStates,
                         int nObservations,
                         int max_obs)
{
    double eps = 1e-30;
    for (int i=0;i<nStates;i++)
        for (int j=0;j<max_obs;j++)
            if (emisP[i][j]==0)
                emisP[i][j]=eps;
    
    for (int i=0;i<nStates;i++)
        for (int j=0;j<nStates;j++)
            if (transP[i][j]==0)
                transP[i][j]=eps;
    
    for (int i=0;i<nStates;i++)
        if (startP[i]==0)
            startP[i]=eps;
    
    double sum;
    for (int i=0;i<nStates;i++)
    {
        sum = 0;
        for (int j=0;j<nStates;j++)
            sum+=transP[i][j];
        for (int j=0;j<nStates;j++)
            transP[i][j]/=sum;
    }
    for (int i=0;i<nStates;i++)
    {
        sum = 0;
        for (int j=0;j<max_obs;j++)
            sum+=emisP[i][j];
        for (int j=0;j<max_obs;j++)
            emisP[i][j]/=sum;
    }
    sum = 0;
    for (int j=0;j<nStates;j++)
        sum+=startP[j];
    for (int j=0;j<nStates;j++)
        startP[j]/=sum;
}//correctModel



static void printModel(vector<vector<double>> &transP,
                       vector<vector<double>> &emisP,
                       vector<double> &startP, int nStates, int nObservations, int max_obs)
{
    std::cout << "\nTRANS: \n";
    for (int r=0;r<nStates;r++)
    {
        for (int c=0;c<nStates;c++)
            std::cout << transP[r][c] << " ";
        std::cout << "\n";
    }
    std::cout << "\nEMIS: \n";
    for (int r=0;r< nStates;r++)
    {
        for (int c=0;c<max_obs;c++)
            std::cout <<"state: "<<r<<","<<c<<" "<<	 emisP[r][c] << " ";
        std::cout << "\n\n";
    }
    std::cout << "\nINIT: \n";
    for (int c=0;c<nStates;c++)
        std::cout << startP[c] << " ";
  
    std::cout << "\n";
}//printModel




//    train(observations,transP, emisP , startP, max_it, nStates,nObservations,max_obs);

/* Calculates maximum likelihood estimates of transition and emission probabilities from a sequence of emissions */
void train(vector<size_t> &observations,
           vector<vector<double>> &transP,
           vector<vector<double>> &emisP,
           vector <double> &startP,
           const int max_iter ,
           int nStates,
           int nObservations,
           int max_obs)
{

    // 1. Initialization
    int iters = 0;
    
    //  correctModel(transP,emisP,startP,     nStates,      nObservations,        max_obs);

    vector<vector<double>> FTRANS(nStates, vector<double>(nStates));
    FTRANS = transP;

    vector<vector<double>> FEMIS(nStates, vector<double>(max_obs  ));
    FEMIS = emisP;

    vector<double> FINIT(nStates);
    FINIT = startP;
    int T= (int) nObservations;

    double logProb = -1 * (std::numeric_limits<double>::max());
    double oldLogProb;
    int data = 0;
    do {
        oldLogProb = logProb;

        // compute a0
        vector<vector<double>> a(nStates, vector<double>(nObservations));
        vector<double> c(nObservations);
        c[0]=0;
        for (int i=0;i<nStates;i++)
        {
            a[i][0] = startP[i]*emisP[i][observations[0]];
            c[0] += a[i][0];
        }
        
        
        // scale the a0(i)
        c[0] = 1/c[0];
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
            startP[i] = YN[i][0];
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
                    if (observations[t] == (size_t) j)
                        numer+= YN[i][t];
                    denom += YN[i][t];
                }
                emisP[i][j] = numer/denom;
            }
        //correctModel(transP,emisP,startP,nStates,nObservations,max_obs);
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
        iters++;
        /* 7. To iterate or not
        data++;
        if (data >= C)
        {
            data = 0;
            iters++;
        }
        */
        
        //FIX ITER AND LOGPROB DIFFERENCE - add R delta
        
    } while (iters<max_iter && logProb>oldLogProb);
   
   // correctModel(FTRANS,FEMIS,FINIT, nStates, nObservations,  max_obs);
    
    printModel(FTRANS,FEMIS,FINIT,(int) nStates, (int) nObservations,(int)max_obs);
    
    /*
    vector<vector<double>> FtransP(nStates, vector<double>(nStates));
    FtransP = FTRANS;

    vector<vector<double>> FemisP(nStates, vector<double>(max_obs  ));
    FemisP = FEMIS;

    vector<double> FstartP(nStates);
    FstartP = FINIT;
     */
    
}//train

//FIX :NEED TO RETURN FxxxP



int main(int argc, const char * argv[])
{
    
    if (argc != 6) {
        std::cerr << "usage: " << argv[0] << " <observation> <Haploid Mean> <output prefix> <max_iter> <delta>" << endl;
        return EXIT_FAILURE;
    }
    
    const double epsi=1e-99;
    const string prefix_file(argv[3]);
    
    const size_t mean=std::stoi(argv[2]);

    //observations--------------------------------------------------------
    vector<size_t> observations;
    const string filename2(argv[1]);
    std::ifstream file;
    size_t inputString;
    file.open(filename2);
    int i=0;
    size_t max=0;
    size_t maxi=0;
    while(file >> inputString){
        observations.push_back(inputString);
        //observations[i]=inputString;
        if(inputString>=max)
            max=inputString;maxi=i;
        i++;
    }
    file.close();
    size_t nObservations=observations.size();
    //----------------------------------------------------------------------
    
    
    
    //max cov value observed or upper cov bound -> max nState---------------
    size_t max_obs = std::min( (int) std::ceil(30.497*mean) , max);
    size_t nStates= (size_t) std::ceil(max_obs/mean);
    
    cout<<"nstates "<<nStates<<endl;
    cout<<"max_cov "<<max_obs<<endl;
    //----------------------------------------------------------------------
 
    
    
    //transition probabilities------------------------------------------------
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
    //---------------------------------------------------------------------------

    
    
    //emission probabilities--------------------------------------------------------
    vector<vector<double>> emisP(nStates, vector<double>(max_obs  ));
    for (size_t i=0;i<nStates;i++)
    {
        for (size_t j=0;j<=max_obs;j++)
        {
            emisP[i][j]=Prpoiss( (int) i ,j, (int) mean );
        }//cout<<endl;
    }
    //---------------------------------------------------------------------------------


    //start probabilities- uniform------------------------------------------------------
    vector<double> startP(nStates);
    for(size_t i=0;i<(nStates);i++)
    {
        startP[i]=1./(nStates);
    }
    //------------------------------------------------------------------------------------
    
    
    
    const int max_it=std::stoi(argv[4]);
    train(observations,transP, emisP , startP, max_it,(int) nStates, (int) nObservations,(int)max_obs);

    std::cout << "\ndone.\n";

    return 0;
}
