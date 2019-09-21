//
//  main.cpp
//  HMM
//
//  Created by RED on 8/26/19.
//  Copyright © 2019 RED. All rights reserved.
//

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <cassert>
#include <cmath>
#include <random>
//#define epsi 1e-99
#include <boost/math/distributions/poisson.hpp> // for normal_distribution
//#include </home/cmb-16/mjc/rdagnew/anaconda3/envs/boost_env/include/boost/math/distributions/fwd.hpp>
using boost::math::poisson;
using boost::math::pdf;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;


/*
static void
read_observation(const string filename, vector<size_t> & T)
{
        std::ifstream file;
//        in(filename);
                     //.c_str());
       // string line;
    //in >> line;
    //float str;
size_t inputString;
    //while ( std::getline(in,str)  ){
//cout<<"read"<<endl;
//    while(in.good()){
//        T.push_back(in.get());
//    }
file.open(filename);
while(file >> inputString)
    T.push_back(inputString);
file.close();

//cout<<T[0]<<endl<<T[1]<<endl<<T[2];
}
*/
double max_over_row(vector<vector<double>> &v , size_t col ,size_t nStates ){

    double maxi=-1 * (std::numeric_limits<double>::max()) ;
//cout<<maxi<<endl;
    //size_t imax=0;
    for(size_t i=0;i< nStates;i++){
        //if ( v[i][col]>maxi ){
            maxi=std::max(v[i][col], maxi);
            //cout<<v[i][col];
            //imax=i;
        //}
    }

//    cout<<" "<<maxi<<endl;//<<v[imax][col]<<endl;
    return maxi;
}
double max_over_rows(vector<vector<double>> &v , size_t col ,vector<vector<double>> &v2 , size_t nextState,size_t nStates ){

    double maxi2=-1 * (std::numeric_limits<double>::max()) ;
//cout<<maxi2<<endl;
    //size_t imax=0;

    //(max(v[, k] + log(hmm$transProbs[, viterbiPath[k + 1] ]    ))

    //max_over_rows(v,k-f,nStates,epsi ) + log(   max_over_row(transP, viterbiPath[(k-f)+1], nStates,epsi)  )

// max_over_rows(v,k-f,transP, viterbiPath[(k-f)+1], nStates,epsi )

    for(size_t i=0;i< nStates;i++){
        //if ( v[i][col]>maxi ){
            maxi2=std::max( v[i][col] + log(v2[i][nextState]) , maxi2);
            //cout<<v[i][col];
            //imax=i;
        //}
    }

    //cout<<" "<<maxi2<<endl;//<<v[imax][col]<<endl;
    return maxi2;
}

double Prpoiss(int cn,  int cov, int Hmean) {
    //Hmean=as.integer(Hmean)
    double result=0;
    if(cn!=0){
        poisson distribution(cn*Hmean);
        result=pdf(distribution, cov);
    }
    else{
        poisson distribution(0.1*Hmean);
        result=pdf(distribution, cov);
    }
    return result;
}

    /*
double Prclip(int cl, int thr ){
    double result = 0;
    if (cl>=thr){
        result=1-eps;
    } else {
        result=eps;
    }
    return result;
}

*/


double emissionPr( size_t index, size_t state, vector<size_t> & observations,size_t mean  ){

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



//vector<float> observations

void viterbi( vector<double> &startP,
    vector<vector<double>> &transP,
    vector<size_t> &observations,
    size_t nStates ,
    size_t mean,
    vector<size_t> &  viterbiPath,
    double epsi ){

    size_t  nObservations  = observations.size();
    vector<vector<double>> v(nStates, vector<double>(nObservations));
    // Init
    for(size_t i=0;i<nStates;i++)
    {
        v[i][0] = log( startP[i] * emissionPr( 1, i , observations,mean) );
    //cout<<"Vi0 "<<v[i][0]<<" emis(1,i,) "<<emissionPr( 1, i , observations,mean)<<endl;
    }
// Iteration
    for(size_t k=1 ; k<nObservations ; k++)
    {
        for(size_t i=0;i<nStates;i++)
        {
            double maxi = -1*epsi;//v[0][k-1] + log(transP[0][i]);
            for(size_t j=0;j<nStates;j++)
            {
                double temp = v[j][k-1] + log(transP[j][i]);
                         //previousState,k-1] + log(hmm$transProbs[previousState,state])
                maxi = std::max(maxi, temp);
                /*
                if (maxi < temp){
                    maxi=temp;
                    imax=
                }
                */
                //cout<<"maxi "<<maxi<<" temp "<<temp<<endl;
            }

            v[i][k] = log(emissionPr( k,i, observations,mean)) + maxi;
            //cout<<"v["<<i<<","<<k<<"]"<<v[i][k]<<endl;

        }
    }
// Traceback
    for(size_t i=0;i<nStates;i++)
    {
        //cout<<max_over_row(v,nObservations-1,nStates,epsi)<<" max row"<<endl;
        if( max_over_row(v,nObservations-1,nStates) == v[i][nObservations-1] )
        {
            viterbiPath[nObservations-1] = i;
    //      cout<<"state n "<<i<<endl;
            break;
        }
    }
    size_t k=nObservations-2;
    for( size_t f=0;f<nObservations-1;   f++ )
    {
        for(size_t i=0;i<nStates;i++)
        {
            //comput max value in column
            //max_over_rows(v,k-f,nStates,epsi ) + log(   max_over_row(transP, viterbiPath[(k-f)+1], nStates,epsi)  )
            double max_rows = max_over_rows(v,k-f,transP, viterbiPath[(k-f)+1], nStates );
            if( max_rows    == v[i][k-f]+ log(transP[i][viterbiPath[(k-f)+1] ] )  )
            {
                viterbiPath[k-f] = i;
                break;
            }
        }
    }


}//viterbi




int main(int argc, const char * argv[]) {

if (argc != 5) {
    std::cerr << "usage: " << argv[0] << " <observation> <nStates> <Haploid Mean> <nObs>" << endl;
    return EXIT_FAILURE;
}

// const string filename(argv[2]);
const string filename2(argv[1]);

// vector<float> observations;

const size_t nStates=std::stoi(argv[2]);

const size_t mean=std::stoi(argv[3]);
const size_t nobs=std::stoi(argv[4]);
//cout<<mean<<nobs<<endl;
vector<size_t> observations(nobs);


//cout<<filename2<<endl;
std::ifstream file;
//        in(filename);
             //.c_str());
// string line;
//in >> line;
//float str;
size_t inputString;
//while ( std::getline(in,str)  ){
//cout<<"read"<<endl;
//    while(in.good()){
//        T.push_back(in.get());
//    }
file.open(filename2);
int i=0;
while(file >> inputString){
    //observations.push_back(inputString);
    observations[i]=inputString;
    i++;
}
file.close();

    //read_observation(filename2, observations);
  //  cout<<filename2<<endl;
    //cout<<observations[0]<<"  "<<observations[1]<<endl;

    size_t nObs=observations.size();

    //cout<<filename2<<" "<<nObs<<endl;

const double epsi=1e-99;
    // start prob
    //uniform(0,nstates-1)
    vector<double> startP(nStates);
    for(size_t i=0;i<nStates;i++){
        startP[i]=0.1;//1/nStates;
//cout<<startP[0];
}

//cout<<"passed"<<endl;
//cout<<epsi<<endl;

    // trans prob
vector<vector<double>> transP(nStates, vector<double>(nStates));
for (size_t i=0;i<nStates;i++){
    for (size_t j=0;j<nStates;j++){
        if(i==j){
            transP[i][j]=1 - ( ( ((double) nStates) -1) * (epsi) );
            cout<<" "<<transP[i][j];
        }
        else{
            transP[i][j]=epsi;
            cout<<" "<<transP[i][j];
        }
    }cout<<endl;
}

//cout<<transP[0][0]<<" "<<transP[1][1]<<" "<<transP[0][1]<<" "<<transP[1][0]<<endl;
//cout<<"before vit";
//viterbi
vector<size_t> viterbiPath(nobs);
viterbi(startP, transP, observations, nStates,mean,viterbiPath, epsi);

std::ofstream output;
output.open("viterout.txt");

for (size_t i=0;i< nObs; i++){
    output<<viterbiPath[i]<<endl;
    //cout<<viterbiPath[i]<<endl;
}
output.close();

//cout<<viterbiPath[0]<<" "<<viterbiPath[nObs-1]<<endl;
return 0;
}//main
