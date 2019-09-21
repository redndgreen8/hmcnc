

#include <iostream>
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


double max_over_row(vector<vector<double>> &v , size_t col ,size_t nStates ){

    double maxi=-1 * (std::numeric_limits<double>::max()) ;
    for(size_t i=0;i< nStates;i++){
            maxi=std::max(v[i][col], maxi);
    }
    return maxi;
}


double max_over_rows(vector<vector<double>> &v , size_t col ,vector<vector<double>> &v2 , size_t nextState,size_t nStates ){
    double maxi2=-1 * (std::numeric_limits<double>::max()) ;

    for(size_t i=0;i< nStates;i++){
            maxi2=std::max( v[i][col] + log(v2[i][nextState]) , maxi2);

    }
    return maxi2;
}

double Prpoiss(int cn,  int cov, int Hmean) {
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




void viterbi( vector<double> &startP,
    vector<vector<double>> &transP,
    vector<size_t> &observations,
    size_t nStates ,
    size_t mean,
    vector<size_t> &  viterbiPath ,size_t nObservations){

    //size_t  nObservations  = observations.size();
    vector<vector<double>> v(nStates, vector<double>(nObservations)  );
    // Init
    for(size_t i=0;i<nStates;i++)
    {
        v[i][0] = log( startP[i] * emissionPr( 1, i , observations,mean) );
    }
// Iteration
    for(size_t k=1 ; k<nObservations ; k++)
    {
        for(size_t i=0;i<nStates;i++)
        {
            double maxi = -1 * (std::numeric_limits<double>::max());
            for(size_t j=0;j<nStates;j++)
            {
                double temp = v[j][k-1] + log(transP[j][i]);
                maxi = std::max(maxi, temp);

            }
            v[i][k] = log(emissionPr( k,i, observations,mean)) + maxi;
        }
    }
// Traceback
    for(size_t i=0;i<nStates;i++)
    {
        if( max_over_row(v,nObservations-1,nStates) == v[i][nObservations-1] )
        {
            viterbiPath[nObservations-1] = i;
            break;
        }
    }
    size_t k=nObservations-2;
    for( size_t f=0;f<nObservations-1;   f++ )
    {
        for(size_t i=0;i<nStates;i++)
        {
            //comput max value in column
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
    std::cerr << "usage: " << argv[0] << " <observation> <nStates> <Haploid Mean> <output prefix>" << endl;
    return EXIT_FAILURE;
}

const string filename2(argv[1]);
const string prefix_file(argv[4]);
const size_t nStates=ceil(std::stof(argv[2]));
//cout<<nStates;
const size_t mean=std::stoi(argv[3]);
//const size_t nobs=std::stoi(argv[4]);
vector<size_t> observations;


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

size_t nObs=observations.size();

const double epsi=1e-99;
vector<double> startP(nStates);

for(size_t i=0;i<nStates;i++){
        startP[i]=1./nStates;
}

    // trans prob
vector<vector<double>> transP(nStates, vector<double>(nStates));
for (size_t i=0;i<nStates;i++){
    for (size_t j=0;j<nStates;j++){
        if(i==j){
            transP[i][j]=1 - ( ( ((double) nStates) -1) * (epsi) );
    //        cout<<" "<<transP[i][j];
        }
        else{
            transP[i][j]=epsi;
    //        cout<<" "<<transP[i][j];
        }
    }//cout<<endl;
}

vector<size_t> viterbiPath(nObs);
viterbi(startP, transP, observations, nStates,mean,viterbiPath,nObs);
string v("viterout.txt");
string filen = prefix_file + "." +v;

std::ofstream output;
output.open(filen.c_str());

for (size_t i=0;i< nObs; i++){
    output<<viterbiPath[i]<<endl;
}
output.close();

return 0;
}//main
