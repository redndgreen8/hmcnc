

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <cassert>
#include <cmath>
#include <boost/math/distributions/poisson.hpp>
using boost::math::poisson;
using boost::math::pdf;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;

static void printModel(vector<vector<double> > &transP,
                       int nStates)
{

    cout<< "\nTRANS: \n";
    for (int r=0;r<nStates;r++)
    {
        cout<<r<<": ";
        for (int c=0;c<nStates;c++)
        {
            cout<< transP[r][c] << " ";
        }
        cout<< "\n";
    }cout<< "\n";
}

double max_over_row(vector<vector<double> > &v , size_t col ,size_t nStates ){

    double maxi=-1 * (std::numeric_limits<double>::max()) ;
    for(size_t i=0;i< nStates;i++){
            maxi=std::max(v[i][col], maxi);
    }
    return maxi;
}


double max_over_rows(vector<vector<double> > &v , size_t col ,vector<vector<double> > &v2 , size_t nextState,size_t nStates ){
    double maxi2=-1 * (std::numeric_limits<double>::max()) ;

    for(size_t i=0;i< nStates;i++){
            maxi2=std::max( v[i][col] + v2[i][nextState] , maxi2);

    }
    return maxi2;
}

double Prpoiss(int cn,  int cov, int Hmean) {
    double result=0;
    const double epsi=1e-99;
    if (Hmean==0){//no alignment in contig
        if(cn!= 0)
            result=epsi;
        else
            result=1-epsi;        
    }
    else if(cov == 10.497*Hmean){//max_obs filtered previously
        if(cn!= 11)
            result=epsi;
        else
            result=1-epsi;
    }
    else if(cn==0){//del_states
        poisson distribution(0.01*Hmean);
        result=pdf(distribution, cov);
    }
    else{
        poisson distribution(cn*Hmean);
        result=pdf(distribution, cov);
    }
    return result;
}




static void correctModel(vector<vector<double> > &transP,
                         int nStates)
{
    double sum;
    for (int i=0;i<nStates;i++)
    {
        sum = 0;
        for (int j=0;j<nStates;j++)
            sum+= std::exp(transP[i][j]);
        for (int j=0;j<nStates;j++)
            transP[i][j]= log(std::exp(transP[i][j])/sum);
    }
}//correctModel


void viterbi( vector<double> &startP,
    vector<vector<double> > &transP,
   vector<vector<double> > &emisP,
    vector<size_t> &observations,
    size_t nStates ,
    size_t mean,
    vector<size_t> &  viterbiPath ,size_t nObservations, size_t max_obs){

    //size_t  nObservations  = observations.size();
    vector<vector<double> > v(nStates, vector<double>(nObservations)  );
    // Init
    size_t obs = std::min(max_obs , observations[0]);
    for(size_t i=0;i<nStates;i++)
    {

        v[i][0] =  startP[i] + emisP[i][obs] ;
    }
// Iteration
    
    for(size_t k=1 ; k<nObservations ; k++)
    {
        size_t obs = std::min(max_obs , observations[k]);
        for(size_t i=0;i<nStates;i++)
        {
            double maxi = -1 * (std::numeric_limits<double>::max());
            for(size_t j=0;j<nStates;j++)
            {
                double temp = v[j][k-1] + transP[j][i];
                maxi = std::max(maxi, temp);

            }
            v[i][k] = emisP[i][obs] + maxi;
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

if (argc != 6) {
    std::cerr << "usage: " << argv[0] << " <coverage observation> <Diploid mean coverage> <output prefix> <scaler> <epsi>" << endl;
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
        i++;
    }
    file.close();

    size_t nObservations=observations.size();

    size_t mean= std::floor( std::stoi(argv[2]) /2 );
    const size_t calc_mean = std::floor( (sum/nObservations)  /2 );

    //max cov value observed or upper cov bound -> max nState---------------
    size_t max_obs = std::min( ( (size_t) std::ceil(10.497*mean)) , max);
    size_t nStates= std::min(  ((int) std::ceil(max_obs/mean)) + 1  ,12   );//+1 zeroth state

    cout<<"nstates "<<nStates<<endl;
    cout<<"max_cov "<<max_obs<<endl;
    cout<<"mean "<<mean<<" calc mean "<<calc_mean<<endl;
    //----------------------------------------------------------------------
double meann;
if (calc_mean==0)
    meann=0.01;
else
    meann = calc_mean;

   // const double epsi=1e-99;
    vector<double> startP(nStates);

    for(size_t i=0;i<(nStates);i++){
            startP[i]=log(1./(nStates));
    }

        // trans prob

        poisson distribution1(3*meann);
        double result3=pdf(distribution1, 3*meann);

        poisson distribution2(2*meann);
        double result2=pdf(distribution2, 3*meann);

        double epsi23 = result2/result3;
/*
        poisson distribution3(mean);
        double result1=pdf(distribution3, mean);

        poisson distribution4(2*mean);
        double result4=pdf(distribution4, mean);

        double epsi21 = result4/result1;

        cout<<"epsi23 "<<epsi23<<" epsi21 "<<epsi21<<endl;
*/

const double epsi = std::pow(0.1,std::stoi(argv[5]));
cout<<"epsi: "<<epsi<<endl;
    //*300
    const double scale = std::stod(argv[4]);
    const double lepsi=log(epsi);

//no epsi
    double beta =  lepsi + ( scale * log(epsi23))  ;
                                    //mean no. of bins for cn=3 call




    vector<vector<double> > transP(nStates, vector<double>(nStates));
    for (size_t i=0;i<nStates;i++)
    {
        for (size_t j=0;j<nStates;j++)
        {
            if(i==j)
            {
                transP[i][j]= log(1 - std::exp(beta) );
            }
            
            
            else
            {
                if ( j==3)
                {
                    transP[i][j]= beta - log(nStates-1) ;

                }
                else
                {
                    transP[i][j]= beta - log(nStates-1);
                
                }
            }
        }
    }

    printModel(transP,nStates);
  //  correctModel(transP,nStates);
   // printModel(transP,nStates);


    vector<vector<double> > emisP(nStates, vector<double>(max_obs+1));
    for (size_t i=0;i<nStates;i++){
        for (size_t j=0;j<=max_obs;j++){
  
                emisP[i][j]=log(Prpoiss( (int) i , j , (int) mean ));

        }
    }



    vector<size_t> viterbiPath(nObservations);


    viterbi(startP, transP,emisP, observations, nStates,mean,viterbiPath,nObservations,max_obs);
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
