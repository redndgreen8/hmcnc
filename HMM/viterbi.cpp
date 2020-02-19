

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


static double max_over_row(vector<vector<double>> &v , size_t col ,size_t nStates ){

    double maxi=-1 * (std::numeric_limits<double>::max()) ;
    for(size_t i=0;i< nStates;i++){
            maxi=std::max(v[i][col], maxi);
    }
    return maxi;
}


static double max_over_rows(vector<vector<double>> &v , size_t col ,vector<vector<double>> &v2 , size_t nextState,size_t nStates ){
    double maxi2=-1 * (std::numeric_limits<double>::max()) ;

    for(size_t i=0;i< nStates;i++){
            maxi2=std::max( v[i][col] + v2[i][nextState] , maxi2);

    }
    return maxi2;
}

double Prpoiss(int cn,  int cov, int Hmean) {
    double result=0;
    const double epsi=1e-99;
    if(cov > 30.497*Hmean){
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

static void correctModel(vector<vector<double>> &transP,
                         int nStates)
{
    double sum;
    for (int i=0;i<nStates;i++)
    {
        sum = 0;
        for (int j=0;j<nStates;j++)
            sum+=transP[i][j];
        for (int j=0;j<nStates;j++)
            transP[i][j]/=sum;
    }
}//correctModel


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
        v[i][0] = log( startP[i] * emissionPr( 0, i , observations,mean) );
    }
// Iteration
    for(size_t k=1 ; k<nObservations ; k++)
    {
        for(size_t i=0;i<nStates;i++)
        {
            double maxi = -1 * (std::numeric_limits<double>::max());
            for(size_t j=0;j<nStates;j++)
            {
                double temp = v[j][k-1] + transP[j][i];
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
            if( max_rows    == v[i][k-f]+ transP[i][viterbiPath[(k-f)+1] ]   )
            {
                viterbiPath[k-f] = i;
                break;
            }
        }
    }


}//viterbi




int main(int argc, const char * argv[]) {

    if (argc < 5) 
    {
        std::cerr << "usage: " << argv[0] << " <coverage observation> <scaler> <output prefix> <calculate_mean/provide_mean[0/1]> <Mean D.coverage>" << endl;
        return EXIT_FAILURE;
    }

    const string prefix_file(argv[3]);
    //const string mean_file(argv[2]);


    //observations--------------------------------------------------------
    vector<size_t> observations;
    const string filename2(argv[1]);
    std::ifstream file;
    size_t inputString;
    file.open(filename2);
    int i=0;
    size_t max=0;
    size_t sum=0;

    //int maxi=0;
    while(file >> inputString)
    {
        observations.push_back(inputString);
        sum=sum+inputString;
        //observations[i]=inputString;
        if(inputString>=max)
        {
            max=inputString;
            //      maxi=i;
        }
    }
    file.close();

    size_t nObservations=observations.size();

    
    size_t mean = std::floor( (sum/nObservations)  /2 );
    if (argv[4]==1){
        mean=std::floor( std::stoi(argv[5]) /2 );
    }



/*
    std::ifstream file2;

    file2.open(mean_file);
    cout<<"to mean"<<endl;
    double inputmean,mean_tmp;
    while(file >> inputmean){
        mean_tmp=inputmean;
        cout<<mean_tmp<<" "<<inputmean<<endl;
    }
    const size_t mean=  std::floor(mean_tmp/2);
    file2.close();

*/

   // size_t nObservations=observations.size();
    //----------------------------------------------------------------------

    //max cov value observed or upper cov bound -> max nState---------------
    size_t max_obs = std::min( ( (size_t) std::ceil(30.497*mean)) , max);
    size_t nStates= ((size_t) std::ceil(max_obs/mean))+1;//+1 zeroth state

    cout<<"nstates "<<nStates<<endl;
    cout<<"max_cov "<<max_obs<<endl;
    cout<<"mean "<<mean<<endl;
    //----------------------------------------------------------------------



    const double epsi=1e-99;
    vector<double> startP(nStates);

    for(size_t i=0;i<(nStates);i++){
            startP[i]=1./(nStates);
    }

        // trans prob

        poisson distribution1(3*mean);
        double result3=pdf(distribution1, 3*mean);

        poisson distribution2(2*mean);
        double result2=pdf(distribution2, 3*mean);

        double epsi23 = result2/result3;

    /*    poisson distribution3(mean);
        double result1=pdf(distribution3, mean);

        poisson distribution4(2*mean);
        double result4=pdf(distribution4, mean);

        double epsi21 = result4/result1;
*/
     //   cout<<"epsi23 "<<epsi23<<" epsi21 "<<epsi21<<endl;


//const double input_epsi = std::stod(argv[5]);

    //*300
    const double scale = std::stod(argv[2]);
    double beta = log(nStates-1) + log(epsi) + (scale * log(epsi23));
                                    //mean no. of bins for cn=3 call





    vector<vector<double>> transP(nStates, vector<double>(nStates));
    for (size_t i=0;i<nStates;i++){
        for (size_t j=0;j<nStates;j++){
            if(i==j)
            {
                transP[i][j]= log(1 - std::exp(beta) );

                //1 - ( ( ((double) (nStates)) -1) * (epsi) );
        //        cout<<" "<<transP[i][j];
            }
            else
            {
                //if (i==2 && j==3)
                transP[i][j]= beta - log(nStates-1);
                //else if(i==2 && j==1)
                //    transP[i][j]= epsi * epsi21;
                ///else
                //    transP[i][j]=epsi;
        //        cout<<" "<<transP[i][j];
            }
        }//cout<<endl;
    }
//    correctModel(transP,nStates);

    vector<size_t> viterbiPath(nObservations);
    viterbi(startP, transP, observations, nStates,mean,viterbiPath,nObservations);
    string v("viterout.txt");
    string filen = prefix_file + "." +v;

    std::ofstream output;
    output.open(filen.c_str());

    for (size_t i=0;i< nObservations; i++){
        output<<viterbiPath[i]<<endl;
    }
    output.close();

return 0;
}//main
