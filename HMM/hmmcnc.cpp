#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include "htslib/faidx.h"
#include <thread>
#include <map>
#include <sstream>
#include <iomanip>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <istream>

using boost::math::poisson;
using boost::math::pdf;
using boost::math::negative_binomial_distribution;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;

int BIN_LENGTH=100;
int MAX_CN=10;
float MISMAP_RATE=0.01;
using namespace std;
double lepsi=-50;

class SNV {
public:
  char refNuc, altNuc;
  int ref;
  int alt;
  int pos;
  SNV(int p, int r, int a, int rc, int ac) : pos(p), refNuc(r), altNuc(a), ref(rc), alt(ac) {}
};
void Reset(vector<vector<vector<double> > > &v) {
  for (int i=0; i < v.size(); i++) {
    for (int j=0; j < v[i].size(); j++) {
      fill(v[i][j].begin(), v[i][j].end(), 0);
    }
  }
}

double PairSumOfLogP(double a, double  b) {
  double res=b;
  if (a!= 0) {
    double m=max(a,b);
    double diff=min(a,b) - m;
    double e=exp(diff);
    double lg=log(1+e);
    res=m + lg;
  }

  return res;
}

double SumOfLogP(vector<double> &vals) {
  if (vals.size() == 0) {
    // Return 0 for error.
    return 0;
  }
  double maxVal = vals[0];
  for (size_t i=1; i < vals.size(); i++) { maxVal=max(maxVal,vals[i]); }
  double expSum=0;
  for (size_t i=0; i < vals.size(); i++) {
    expSum+= exp(vals[i]-maxVal);
  }
  return maxVal+log(expSum);
}


using namespace std;
class ThreadInfo {
public:
  htsFile *htsfp;
  hts_idx_t *bamidx;
  bam_hdr_t *samHeader;			
  faidx_t *fai;
  int *lastSeq;
  string refFileName;
  pthread_mutex_t *semaphore;
  vector<string> *contigNames;
  vector<int>    *contigLengths;
  vector<int> procChroms;
  vector<vector<int> > *covBins;
  vector<vector<SNV> > *snvs;
  vector<vector<int> > *copyNumber;
  vector<vector<double> > *transP, *emisP;
  vector<double> *startP;
  int maxCov, maxState;
  double mean;
  double var;
  double lepsi;
  double scale;
};

//
// Not super 
void ReadCoverage(string &covFileName,
		  vector<string> &contigNames,
		  vector<vector<int> > &covBins) {
  ifstream covFile(covFileName.c_str());
  string chrom="", curChrom;
  int start, end;
  int cov;
  int last=0;
  int length;
  covFile.seekg(0, std::ios::end);    // go to the end
  length = covFile.tellg();           // report location (this is the length)
  covFile.seekg(0, std::ios::beg);    // go back to the beginning
  char *buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  covFile.read(buffer, length);       // read the whole file into the buffer
  cerr << "read cov buffer of len " << length << endl;
  covFile.close();  
  int i=0;
  string contigName("");
  int curContig=0;
  if (length > 0) {
    covBins.push_back(vector<int>() );
  }
  while (i < length) {
    while (i < length and isspace(buffer[i])) { i++; }
    int c=i;    
    while (i < length and isspace(buffer[i]) == false)  { i++; }
    if (i < length) {
      if (i-c > contigName.size()) { contigName.resize(i-c);}
      sscanf(&buffer[c], "%s", contigName.c_str());
      sscanf(&buffer[i], "%d	%d	%d", &start, &end, &cov);
      if (contigName != contigNames[curContig]) {
	covBins.push_back(vector<int>());
	curContig++;	
      }
      covBins[curContig].push_back(cov);
    }
    while (i < length and buffer[i] != '\n') { i++;};
  }
}

void WriteCovBed(string &covFileName,
		 vector<string> &contigNames,
		 vector<vector<int> > &covBins) {
  ofstream covFile(covFileName.c_str());
  for (int c=0; c < contigNames.size(); c++) {
    for (int i=0; i < covBins[c].size(); i++) {
      covFile << contigNames[c] << "\t" << i*100 << "\t" << (i+1)*100 << "\t" << covBins[c][i] << endl;
    }
  }
}

static void printModel(vector<vector<double> > &transP)
{
  ostream* strm =&cout;
  *strm << "\nTRANS: \n";
  for (int r=0;r<transP.size();r++)
    {
      *strm << r <<":";
      for (int c=0;c<transP[r].size();c++)
        {
	  *strm << "\t"<< transP[r][c];
        }
      *strm << endl;
    }
  *strm << endl;
}

static void printEmissions(vector<vector<double> > &emisP)

{
  cerr << "Emissions matrix: \n";
  for (size_t r=0;r<emisP.size(); r++) 
    {
      for (size_t c=0;c<emisP[r].size(); c++ )
	{
	  cerr << std::setw(4) << emisP[r][c];
	  if (c+1 < emisP[r].size()) 
	    {
	      cerr << " ";
	    }	  
	}
      cerr << endl;
    }
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

double LgNegBinom(int cn, int cov, float Hmean, float Hvar) {
  double result=0;
  float r, p;

  if (Hmean == 0 or Hvar == 0) {
    return 0;
  }
  //
  // Since actual variance is unknown at high copy number states, assume linear increase.
  //
      
  if (Hmean==0) {//no alignment in contig
    if (cn!= 0)
      result=lepsi;
    else
      result=0;    
  }
  else if (cov >= MAX_CN*Hmean) {//max_obs filtered previously
    if(cn!= MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    //
    // Use poisson for 0-state assuming it's a random mismap process.
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else 
      {
	result=log(prob);
      }
  }
  else{
    Hmean*=cn;
    Hvar*=cn;
    p=Hmean/Hvar;
    r=Hmean*(Hmean/Hvar)/(1-Hmean/Hvar);
    
    negative_binomial_distribution<double> distribution(r,p);    
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else
      {
	result=log(prob);

      }
  }
  return result;  
}

double LgPrpoiss(int cn,  int cov, int Hmean) {
  double result=0;

  // Boundary case
  if (Hmean==0){//no alignment in contig
    if(cn!= 0)
      // Cannot emit 0 reads from non-zero states.
      result=lepsi;
    else
      // Can only emit 0
      result=0;
  }
  else if (cov >= MAX_CN*Hmean){//max_obs filtered previously
    //
    // Have state that catches ultra-high probability states.
    if(cn != MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) {
      result = lepsi;      
    }
    else {
      result = log(prob);
    }
  }
  else {
    poisson distribution(cn*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0)
      {
	result=lepsi;
      }
    else 
      {
	result=log(prob);	
      }
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

void ForwardBackwards( vector<double> &startP,
		       vector<vector<double> > &transP,
		       vector<vector<double> > &emisP,
		       vector<int> &obs,
		       vector<vector<double> > &f,
		       vector<vector<double> > &b) {
  int nObs=obs.size();
  int nStates=startP.size();
  //
  // Initialize first col from standing.
  //

  // f and b are of length nObs+1
  // The first/last prob from an observation are at:
  //  f[1], f[nObs], from obs[0...nObs-1]
  //  and b[nObs], b[1] from obs[nObs-1 ... 0]
  // 
 
  for (int j=0; j < nStates; j++) {
    f[j][0] = startP[j];
  }
  vector<double> col(nStates, 0);
  //
  // Shouldn't have 0 states. Even if the coverage is empty for
  // an entire chrom, should have 0 state.
  //
  assert(nStates > 0);  
  //
  // If just one state (e.g. all zeros), same prob.
  //

  if (nStates == 1) {
    fill(f[0].begin(), f[0].end(), -1);
    fill(b[0].begin(), b[0].end(), -1);
    return;
  }
  
  for (int k=0; k < nObs; k++) {
    cout << "f: " << k << "\t" << obs[k];
    for (int i=0; i < nStates; i++) {
      double colSum=0;
      for (int j=0; j < nStates; j++) {
	colSum = PairSumOfLogP(colSum, f[j][k] + transP[j][i]);
      }
      f[i][k+1] = colSum + emisP[i][obs[k]];
      cout << "\t" << f[i][k+1];
    }
    cout << endl;
  }

  // back
  for (int j=0; j < nStates; j++) {
    b[j][nObs] = startP[j];
  }

  for (int k=nObs; k > 0; k--) {
    for (int i=0; i < nStates; i++) {
      double colSum=0;
      for (int j=0; j < nStates; j++) {
	colSum = PairSumOfLogP(colSum, b[j][k] + transP[j][i]);
      }
      b[i][k-1] = colSum + emisP[i][obs[k-1]];
    }      
  }

  for (int k=0; k< nObs; k++) {
    double max=0;
    double sum=0;
    int maxi=0;
    double pSum=0;
    cout << k << "\t" << obs[k];
    for (int s=0; s < nStates; s++) {
      pSum = f[s][k] + b[s][k+1];
      cout << "\t" << pSum;
      sum = PairSumOfLogP(sum, pSum);
      if (max == 0 or max < pSum) {
	max=pSum;
	maxi=s;
      }
    }

    if (sum !=  0) {
      cout << "\tPD:\t" << maxi << "\t" << exp(max-sum) << "\t" << obs[k];
    }
    cout << endl;
  }
}

void BaumWelchEMStep(vector<double> &startP,
		     vector<vector<double> > &transP,
		     vector<vector<double> > &emisP,
		     vector<int> &obs,
		     vector<vector<double > > &f,
		     vector<vector<double> > &b,
		     vector<vector<double> > &updateTransP,
		     vector<vector<double> > &updateEmisP) {
  int nStates = startP.size();
  ForwardBackwards( startP, transP, emisP, obs, f, b);

  vector<vector<double > > expTransP, expEmisP;
  expTransP.resize(transP.size());
  expEmisP.resize(emisP.size());


  //
  // E step
  //
  for (int j=0; j < transP.size(); j++) {
    expTransP[j].resize(transP[j].size());
    for (int k=0; k < transP[j].size(); k++) {
      double logSum=0;
      for (int i=0; i< obs.size()-1; i++) {
	assert(j < f.size());
	assert(i < f[j].size());
	assert(j < transP.size());
	assert(k < transP[j].size());
	assert(j < emisP.size());
	assert(obs[i] < emisP[j].size());
	//	cout << "BW: " << j << " " << k << "\t" << logSum << "\ttot: " << f[j][i] + transP[j][k] + emisP[j][obs[i+1]] + b[k][i+1] << "\t" << "\tf: " << f[j][i] << "\tt: " << transP[j][k] << "\te: " << emisP[j][obs[i+1]] << "\tb: " << b[k][i+1] << "\t";
	logSum = PairSumOfLogP(logSum, f[j][i] + transP[j][k] + emisP[j][obs[i+1]] + b[k][i+1]);
      }
      expTransP[j][k] = logSum;
    }
  }
  for (int k=0; k < nStates; k++) {
    expEmisP[k].resize(emisP[k].size());
    fill(expEmisP[k].begin(), expEmisP[k].end(), 0);
    for (int i=0; i < obs.size(); i++) {
      expEmisP[k][obs[i]] = PairSumOfLogP(expEmisP[k][obs[i]], f[k][i] + b[k][i+1]);
    }    
  }
  //
  // M step.
  //
  updateTransP.resize(nStates);
  for (int j=0; j < nStates; j++) {
    double colSum=0;
    updateTransP[j].resize(nStates);
    for (int k=0; k< nStates; k++) {
      colSum=PairSumOfLogP(colSum, expTransP[j][k]);
    }
    for (int k=0; k< nStates; k++) {
      updateTransP[j][k] = expTransP[j][k] - colSum;
    }
  }
  updateEmisP.resize(nStates);
  for (int j=0; j < nStates; j++) {
    updateEmisP[j].resize(emisP[j].size());
    double denom=0;
    for (int i=0; i < expEmisP[j].size(); i++) {
      denom = PairSumOfLogP(denom, expEmisP[j][i]);
    }

    for (int i=0; i < expEmisP[j].size(); i++) {
      if (denom > 0) {
	updateEmisP[j][i] = expEmisP[j][i] - denom;
      }
      else {
	updateEmisP[j][i] = 0;
      }
    }
  }  
}
		
  

void viterbi( vector<double> &startP,
	      vector<vector<double> > &transP,
	      vector<vector<double> > &emisP,
	      vector<int> &observations,
	      size_t mean,
	      vector<int> & viterbiPath, int maxAllowedCov){

  //size_t  nObservations  = observations.size();
  size_t nObservations=observations.size();
  int nStates=transP.size();
  vector<vector<double> > v(nStates, vector<double>(observations.size()) );
  //  vector<vector<double> > f(nStates, vector<double>(observations.size()) );
  //  vector<vector<double> > b(nStates, vector<double>(observations.size()) );
  vector<vector<double> > opt(nStates, vector<double>(observations.size())  );  
  // Init
  int obs = std::min(maxAllowedCov , observations[0]);
  int last=observations.size();
  for(size_t i=0;i<nStates;i++)
    {
      v[i][0] = startP[i] + emisP[i][obs] ;
      //      f[i][0] = startP[i] + emisP[i][obs] ;
      //      b[last-i-1][0] = startP[i] + emisP[i][obs] ;
    }
  // Iteration
    
  for(size_t k=1 ; k<nObservations ; k++)
    {
      size_t obs = std::min(maxAllowedCov, observations[k]); 
      for(size_t i=0;i<nStates;i++)
        {
	  double maxProb = v[0][k-1] + transP[0][i];
	  int maxState=0;
	  for(size_t j=1;j<nStates;j++)
            {
	      double rowProb = v[j][k-1] + transP[j][i];
	      if (rowProb > maxProb) {
		maxState=j;
		maxProb=rowProb;
	      }
            }
	  v[i][k] = emisP[i][obs] + maxProb;
	  opt[i][k] = maxState;	  
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
  size_t lastObservation=nObservations-2;
  //
  // Find highest-scoring final entry.
  //
  size_t rowIndex=nObservations-2;
  double maxScore=v[0][rowIndex];
  int    maxState=0;
  for (size_t i=1;i < nStates; i++) {
    if (maxScore < v[i][rowIndex]) {
      maxState=i;
      maxScore=v[i][rowIndex];
    }
  }
  viterbiPath[lastObservation] = maxState;
  
  
  for( size_t f=lastObservation; f > 0; f--) {
      viterbiPath[f] = maxState;
      maxState=opt[maxState][f];
    }

}//viterbi

class CountNuc {
public:
  int index;
  int count;
  int operator<(const CountNuc &rhs) {
    return count < rhs.count;
  }
};

int StoreSNVs(char *contigSeq,
	      int contigLength,
	      float mean,
	      vector<int> &nA, vector<int> &nC, vector<int> &nG, vector<int> &nT, vector<int> &nDel,
	      vector<SNV> &snvs) {
  vector<int* > fPtr(5);
  //
  // Easier for typing
  //
  const char *nucs="ACGTd";
  fPtr[0] = &nA[0];
  fPtr[1] = &nC[0];
  fPtr[2] = &nG[0];
  fPtr[3] = &nT[0];
  fPtr[4] = &nDel[0];
  vector<CountNuc > counts;
  counts.resize(5);
  for (int i =0; i < contigLength; i++) {
    for (int n=0; n < 5; n++) {
      counts[n].index=n;
      counts[n].count = fPtr[n][i];
    }
    sort(counts.begin(), counts.end());
    char refNuc=toupper(contigSeq[i]);
    if (counts[4].index != 4 and
	counts[3].index != 4 and
	refNuc != 'N' and
	counts[3].count > 0.25*mean and
	counts[4].count > 0.25*mean )  {
      //
      // Top most frequent are not a deletion.
      //
      if (nucs[counts[4].index] == refNuc) {
	snvs.push_back(SNV(i, refNuc, nucs[counts[3].index], counts[4].count, counts[3].count));
      }
      else if (nucs[counts[3].index] == refNuc) {
	snvs.push_back(SNV(i, refNuc, nucs[counts[4].index], counts[4].count, counts[3].count));
      }

    }
  }
  return 1;
}


int IncrementCounts(bam1_t *b,
		    int contigLength,
		     vector<int> &nA, vector<int> &nC, vector<int> &nG, vector<int> &nT, vector<int> &nDel) {
  int readLength = b->core.l_qseq;			
  if (readLength < BIN_LENGTH or b->core.qual < 10 or b->core.flag & 0x800) {
    return 0;
  }
	  
  vector<char> seq(readLength);
  uint8_t *q = bam_get_seq(b);
  for (int i=0; i < readLength; i++) {seq[i]=seq_nt16_str[bam_seqi(q,i)];	}
  uint32_t* cigar = bam_get_cigar(b);
  int refLen = bam_cigar2rlen(b->core.n_cigar, cigar);

  int qPos=0;
  int refPos = b->core.pos;
  int ci;
  int regionOffset=refPos;
  bool first=true;
  for (ci=0; ci < b->core.n_cigar; ci++) {
    int opLen=bam_cigar_oplen(cigar[ci]);
    int op=bam_cigar_op(cigar[ci]);
	  
    if (op == BAM_CSOFT_CLIP) {
      qPos += opLen;
      continue;
    }
    if (op == BAM_CINS) {
      qPos += opLen;
      continue;
    }
    if (op == BAM_CDEL) {
      int stop=refPos+opLen;
      for (; refPos < stop and refPos < contigLength; refPos++) {
	nDel[regionOffset]+=1;
	regionOffset++;
      }
      continue;
    }
    if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) {
      if (refPos + opLen <= 0) {
	refPos += opLen;
	qPos += opLen;
	continue;
      }
      else {
	for (int p=0; p < opLen; p++) {
	  if (refPos >= contigLength){ 
	    break;
	  }
	  if (refPos >= 1) {
	    first=false;
	    char nuc=toupper(seq[qPos]);
		    
	    assert(regionOffset < nA.size());
	    if (nuc == 'A') { nA[regionOffset]++;}
	    if (nuc == 'C') { nC[regionOffset]++;}
	    if (nuc == 'G') { nG[regionOffset]++;}
	    if (nuc == 'T') { nT[regionOffset]++;}
	    regionOffset++;
	  }
	  refPos++;
	  qPos++;
	}
      }
    }
  }
  return 1;
}


void ParseChrom(ThreadInfo *threadInfo) {

  while (*(threadInfo->lastSeq) < (*(*threadInfo).contigNames).size()) {
    //
    // Grab current chrom to process
    //
    pthread_mutex_lock(threadInfo->semaphore);

    int curSeq = *((*threadInfo).lastSeq);
    *(threadInfo->lastSeq) = *(threadInfo->lastSeq) + 1;


    //
    // Deal with race condition by double checking curSeq;
    //
    if (curSeq >= threadInfo->contigNames->size()) {
      pthread_mutex_unlock(threadInfo->semaphore);      
      break;
    }

    (*threadInfo).procChroms.push_back(curSeq);
    (*threadInfo).covBins->push_back(vector<int>());
    (*threadInfo).snvs->push_back(vector<SNV>());
    pthread_mutex_unlock(threadInfo->semaphore);
    
    int contigLength=threadInfo->contigLengths->at(curSeq);


    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);
    
    stringstream regionStrm;
    regionStrm << (*(*threadInfo).contigNames)[curSeq];// << ":1-" << contigLength;

    string region=regionStrm.str();


    hts_itr_t *regionIter = sam_itr_querys(threadInfo->bamidx, threadInfo->samHeader, region.c_str());
    int chromLen;
    char *chromSeq = fai_fetch(threadInfo->fai, region.c_str(), &chromLen);
    
    bool continueParsing=true;
    vector<bam1_t*> reads; //(bam_init1());
    long totalSize=0;
    int chunkNumber=0;
    int totalReads=0;
    int endpos=0;
    while (continueParsing) {
      int totalSize=0;
      reads.resize(0);
      pthread_mutex_lock(threadInfo->semaphore);
      int bufSize=0;
      while (bufSize < 500000000 and continueParsing) {
	bam1_t *b = bam_init1();	
	int res=sam_itr_next(threadInfo->htsfp, regionIter, b);
	bufSize+= b->l_data;
	totalSize+= b->l_data;
	
	if (res < 0) { // or totalReads < 15000) {
	  continueParsing = false;
	  cerr << "Ending parsing of " << region << " with " << totalSize << " data and " << chunkNumber << " iterations on res " << res << endl;
	  break;
	}
	endpos=bam_endpos(b);
	reads.push_back(b);
	++totalReads;
      }
      cerr << "Reading contig index " << curSeq << " (chunk " << chunkNumber << ")\t" << (*threadInfo->contigNames)[curSeq] << " " << reads.size() << "\t" << totalReads << " ending at " << endpos << endl;
      ++chunkNumber;
      pthread_mutex_unlock(threadInfo->semaphore);

      for (int readIndex=0; readIndex < reads.size(); readIndex++) {
	bam1_t *b=reads[readIndex];
	IncrementCounts(b, contigLength, nA, nC, nG, nT, nDel);
	bam_destroy1(b);
      }
    }
    // Never compute in the last bin
    int nBins=contigLength/BIN_LENGTH;

    for (int bin=0; bin < nBins; bin++) {
      int start=bin*BIN_LENGTH;
      int end=min((bin+1)*BIN_LENGTH, contigLength);
      int totCov=0;
      for (int bp=start; bp < end; bp++) {
	totCov+=nA[bp] + nC[bp] + nG[bp] + nT[bp] + nDel[bp];
      }
      (*threadInfo->covBins)[curSeq][bin] =totCov/BIN_LENGTH;
    }
    
    //
    // Detect SNVS
    //        
    StoreSNVs(chromSeq, chromLen, threadInfo->mean,
	      nA, nC, nG, nT, nDel,
	      (*threadInfo->snvs)[curSeq]);
    cerr << "Stored " << (*threadInfo->snvs)[curSeq].size() << " snvs " << endl;
    
  }
  pthread_exit(NULL);    
}
    
/*
    if ((*threadInfo->covBins)[curSeq].size() > 0) { 
      //
      double chromMean;
      long chromTot=0;
      viterbi( *threadInfo->startP, *threadInfo->transP, *threadInfo->emisP, (*threadInfo->covBins)[curSeq], threadInfo->mean, (*threadInfo->copyNumber)[curSeq], threadInfo->maxCov);      
    }
*/

vector<int> NucMap;
int GetRefAndAlt(char refNuc, vector<int> &counts, int &ref, int &alt) {
  vector<int> sCounts=counts;
  sort(sCounts.begin(), sCounts.end());
  
  int second = sCounts[2];
  int first  = sCounts[3];
  int firstIndex=-1, secondIndex=-1;
  for (int i=0; i < 4; i++) {
    if (firstIndex == -1 and counts[i] == first) {
      firstIndex=i;
    }
    else if (secondIndex == -1 and counts[i] == second) {
      secondIndex=i;
    }
  }
  if (firstIndex == -1 or secondIndex == -1) {
    ref=0; alt=0;
    return 4;
  }
  int refNucIndex=NucMap[refNuc];
  if (first == 0 or second/first < 0.2) {
    ref=0;
    alt=0;
    return 4;
  }
  if (firstIndex == refNucIndex) {
    ref=first;
    alt=second;
    return secondIndex;
  }
  else {
    ref=second;
    alt=first;
    return firstIndex;
  }  
}

const char* nucs = "ACGTN";

static int pileup_blank(void *data, bam1_t *b) {
  return 0;
}
int EstimateCoverage(string &bamFileName, vector<vector<int> > &allCovBins, vector<string> &chroms, vector<int> &lengths, string &useChrom, double &mean, double &var) {
  int useChromIndex=0;
  if (useChrom == "") {
    int maxLen=0;
    for (int i=0; i < lengths.size(); i++) {
      if (lengths[i] > maxLen) {
	useChrom = chroms[i];
	maxLen=lengths[i];
	useChromIndex=i;
      }
    }
  }
  cerr << "Estimating coverage from " << useChrom << endl;
  int contigLength=0;
  for (int i=0; i < chroms.size(); i++) {
    if (chroms[i] == useChrom) {
      contigLength=lengths[i];
      useChromIndex=i;
      break;
    }
  }
  
  if (contigLength == 0) {
    cerr << "ERROR Could not estimate coverage." <<endl;
    exit(1);
  }
  if (allCovBins.size() > 0) {
    assert(allCovBins[useChromIndex].size() == lengths[useChromIndex]/100);
    int lastBin=allCovBins[useChromIndex].size();
    if (lastBin == 0) {
      cerr << "ERROR. Could not estimate coverage using precomputed bins." << endl;
      exit(1);
    }
    long totCov=0;
    for (int binIndex=0; binIndex<lastBin; binIndex++) {
      totCov+=allCovBins[useChromIndex][binIndex];
    }
    mean=totCov/lastBin;
    //
    // Recompute summary stats using limited data
    int nSamples=0;
    totCov=0;
    long totCovSq=0;
    for (int binIndex=0; binIndex<lastBin; binIndex++) {
      if (allCovBins[useChromIndex][binIndex] > 0.25*mean and allCovBins[useChromIndex][binIndex] < 1.75 * mean) {
	totCov+=allCovBins[useChromIndex][binIndex];
	totCovSq+=allCovBins[useChromIndex][binIndex]*allCovBins[useChromIndex][binIndex];
	nSamples++;
      }
    }
    if (nSamples > 0) {
      mean=totCov/nSamples;
      var=totCovSq/nSamples-mean*mean;
      cerr << "Estimating coverage on precomputed bins " << nSamples << " ending at " << mean << "\t" << var << endl;
    }
    else {
      cerr << "Could not estimate coverage using precomputed coverage on chrom " << useChrom << endl;
      exit(1);
    }
  }
  else {
    htsFile *htsfp;
    vector<int> covBins;
    htsfp = hts_open(bamFileName.c_str(),"r");
  
    hts_idx_t *bamidx;
    if ((bamidx = sam_index_load(htsfp, bamFileName.c_str())) == 0) {
      cerr << "ERROR reading index" << endl;
      exit(0);
    }
    
    const htsFormat *fmt = hts_get_format(htsfp);
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cout << "Cannot determine format of input reads." << endl;
      exit(1);
    }

    bam_hdr_t *samHeader;			
    samHeader = sam_hdr_read(htsfp);
    
    hts_itr_t *regionIter = sam_itr_querys(bamidx, samHeader, useChrom.c_str());

    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);
  
  
    bool continueParsing=true;
    int nSamples=0;
    
    int curEndPos=0;
    int curCovBin=0;
    long totalSize;
    
    while (continueParsing) {    
      int bufSize=0;
      int nReads=0;
      while (bufSize < 100000000 and continueParsing) {
	bam1_t *b = bam_init1();	
	int res=sam_itr_next(htsfp, regionIter, b);
	bufSize+= b->l_data;
	totalSize+= b->l_data;
	
	if (res < 0) {
	  continueParsing = false;
	  break;
	}
	if (IncrementCounts(b, contigLength, nA, nC, nG, nT, nDel)) {
	  curEndPos=bam_endpos(b);
	}
	bam_destroy1(b);
	++nReads;
      }

      //
      // Compute coverage for bins
      int lastBin=(curEndPos-30000)/100;
      for (int binIndex=curCovBin; binIndex < lastBin; binIndex++) {
	int binTot=0;
	for (int nuc=binIndex*100; nuc < (binIndex+1)*100; nuc++) {
	  binTot += nA[nuc] + nC[nuc]+ nG[nuc] + nT[nuc] + nDel[nuc];
	}
	covBins.push_back(binTot/100);
	
      }
      curCovBin=lastBin;
    //
    // Get summary statistics
    //
      double totCov=0;
      double totCovSq=0;
      //
      // First pass gets close to CN=2
      //
      
      for (int binIndex=0; binIndex<lastBin; binIndex++) {
	totCov+=covBins[binIndex];
      }

      mean=totCov/lastBin;
      totCov=0;
      nSamples=0;
      for (int binIndex=0; binIndex<lastBin; binIndex++) {
	if (covBins[binIndex] > 0.25*mean and covBins[binIndex] < 1.75 * mean) {
	  totCov+=covBins[binIndex];
	  totCovSq+=covBins[binIndex]*covBins[binIndex];
	  nSamples++;
	}
      }
      if (nSamples > 0) {
	mean=totCov/nSamples;
	var=totCovSq/nSamples-mean*mean;
	cerr << "Estimating coverage " << nReads << " ending at " << curEndPos << "\t" << mean << "\t" << var << endl;
      }
      if (nSamples > 80000) {
	return 1;
      }
    }
  }
  return 1;
}






int main(int argc, const char* argv[]) {
  int nproc=4;
  double scale=2;
  
  if (argc < 3) {
    cout << "usage: hmmcnc reference.fa" << endl
         << "   -a alignments    Read alignments from this file and calculate depth on the fly." << endl
         << "   -b bed           Read depth bed from this file. Skip calculation of depth." << endl
         << "   -S snv-file      Read SNVs from this file (when not estimating from a BAM)" << endl
	 << " Options controlling depth calculation " << endl      
	 << "   -e value (float)   Value of log-epsilon (-500)." << endl      
	 << "   -S value (float)   Scalar for transition probabilities (pre Baum-Welch) (10)" << endl
	 << "   -m value [pois|nb] Coverage model to use, Poisson (pois), or negative binomial (nb). Default nb." << endl
	 << "   -x value Max state to allow (10)" << endl
	 << " -t value (int)     Number of threads (4) " << endl            
	 << " -c contig          Use this contig to estimate coverage. By default, longest contig." << endl      
	 << " Options controlling output:" << endl
         << " -o file            Output to this file (stdout)." << endl
	 << " -C contig          Only run hmm on this chrom." << endl
         << " -B bed             Write coverage bed to this file." << endl
	 << " -M (flag)          Merge consecutive bins with the same copy number." << endl;
    exit(1);
  }
  NucMap.resize(256,4);
  NucMap[(int)'A']=0;
  NucMap[(int)'C']=1;
  NucMap[(int)'G']=2;
  NucMap[(int)'T']=3;	
  int maxState=10;	
  string bamFileName="";
  string referenceName=argv[1];
  typedef enum { POIS, NEG_BINOM  } MODEL_TYPE;
  MODEL_TYPE model=NEG_BINOM;
  string useChrom="";
  string hmmChrom="";
  string covBedInFileName="", covBedOutFileName="";
  bool   mergeBins=false;
  string outBedName="";
  string outFileName="";
  string snvFile="";
  if (argc > 2) {
    int argi=2;
    while (argi < argc) {
      if (strcmp(argv[argi], "-a") == 0) {
	++argi;
	bamFileName=argv[argi];
      }
      else if (strcmp(argv[argi], "-s") == 0) {
	++argi;
	snvFile=argv[argi];
      }      
      else if (strcmp(argv[argi], "-t") == 0) {
	++argi;
	nproc=atoi(argv[argi]);
      }
      else if (strcmp(argv[argi], "-b") == 0) {
	++argi;
	covBedInFileName = argv[argi];
      }
      else if (strcmp(argv[argi], "-B") == 0) {
	++argi;
	covBedOutFileName = argv[argi];
      }      
      else if (strcmp(argv[argi], "-e") == 0) {
	++argi;
	lepsi=atof(argv[argi]);
      }
      else if (strcmp(argv[argi], "-s") == 0) {
	++argi;
	scale=atof(argv[argi]);
      }
      else if (strcmp(argv[argi], "-m") == 0) {
	++argi;
	if (strcmp(argv[argi], "pois") == 0) {
	  model=POIS;
	}
      }
      else if (strcmp(argv[argi], "-o") == 0) {
	++argi;
	outFileName = argv[argi];
      }      
      else if (strcmp(argv[argi], "-c") == 0) {
	++argi;
	useChrom = argv[argi];
      }
      else if (strcmp(argv[argi], "-C") == 0) {
	++argi;
	hmmChrom = argv[argi];
      }
      else  if (strcmp(argv[argi], "-M") == 0) {
	mergeBins=true;
      }
      else if (strcmp(argv[argi], "--earlyExit") == 0) {
	++argi;
      }
      ++argi;
    }
  }

  string faiFileName=referenceName + ".fai";
  ifstream faiIn(faiFileName.c_str());
  if (faiIn.good() == false) {
    cerr << "ERROR. Reference is not indexed, or could not open .fai file" << endl;
    exit(1);
  }

  vector<string> contigNames, allContigNames;
  vector<int>    contigLengths, allContigLengths;
  while (faiIn) {
    string line;
    getline(faiIn, line);
    stringstream strm(line);
    if (line != "") {
      string contig;
      int length;
      
      strm >> contig;
      strm >> length;      
      allContigNames.push_back(contig);
      allContigLengths.push_back(length);
      if (hmmChrom == "" or contig == hmmChrom) {
	contigNames.push_back(contig);
	contigLengths.push_back(length);
      }
    }
  }
  vector<vector<int> > covBins;
  double mean;
  double var;
  if (covBedInFileName != "") {
    ReadCoverage(covBedInFileName, contigNames, covBins);
  }

  
  EstimateCoverage(bamFileName, covBins, allContigNames, allContigLengths, useChrom, mean, var);
  
  cerr << "Got cov " << mean << "\t" << var << endl;
  
  //
  // Get the header of the bam file
  //
  htsFile *htsfp=NULL;
  bam_hdr_t *samHeader=NULL;
  hts_idx_t *bamidx=NULL;

  
  if (bamFileName != "") {
    htsfp = hts_open(bamFileName.c_str(),"r");
    samHeader = sam_hdr_read(htsfp);
    if ((bamidx = sam_index_load(htsfp, bamFileName.c_str())) == 0) {
      cerr << "ERROR reading index" << endl;
      exit(0);
    }
    const htsFormat *fmt = hts_get_format(htsfp);
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cout << "Cannot determine format of input reads." << endl;
      exit(1);
    }    
  }

  faidx_t *fai = fai_load_format(referenceName.c_str(), FAI_FASTA);

  // 
  // Read index for random io
  //

  pthread_t *threads = new pthread_t[nproc];
  vector<ThreadInfo> threadInfo(nproc);		
  pthread_mutex_t semaphore;		
  pthread_mutex_init(&semaphore, NULL);
  pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
  int curSeq=0;

  vector<vector<SNV> > snvs;
  vector<vector<int > > copyNumber;
  vector<vector< vector<double> > > f,b;
  copyNumber.resize(contigLengths.size());

  if (covBedInFileName == "") {
    covBins.resize(contigLengths.size());
    for (int c=0; c < contigLengths.size(); c++ ) {
      covBins[c].resize(contigLengths[c]/BIN_LENGTH);
      copyNumber[c].resize(contigLengths[c]/BIN_LENGTH);       
    }
  }


  //max cov value observed or upper cov bound -> max nState---------------
  size_t nStates= std::min( maxState , MAX_CN  ) + 1; //+1 zeroth state
  MAX_CN=nStates+1;
  cerr << "Computing copy-number up to " << MAX_CN << endl;
  //----------------------------------------------------------------------

  if (mean == 0) {
    std::cerr << "mean is zero, Exiting" << endl;
    return EXIT_FAILURE;
  }

  vector<double> startP(nStates);

  for(size_t i=0;i<(nStates);i++){
    startP[i]=log(1./(nStates));
  }

  f.resize(covBins.size());
  b.resize(covBins.size());
  for (size_t c=0; c < covBins.size(); c++) {
    f[c].resize(nStates);
    b[c].resize(nStates);
    for (size_t s=0; s < nStates; s++) {
      f[c][s].resize(covBins[c].size() + 1);
      b[c][s].resize(covBins[c].size() + 1);
    }
  }
  
  // trans prob, scale 3->2 by overlap of pdf.

  poisson distribution1(3*mean/2);
  double result3=pdf(distribution1, 3*mean/2);

  poisson distribution2(2*mean/2);
  double result2=pdf(distribution2, 3*mean/2);

  double epsi23 = result2/result3;

  //*300

  //no epsi
  double small=-2;
  //  double beta =  small + ( scale * log(epsi23))  ;
  double beta=small;
  //mean no. of bins for cn=3 call

  vector<vector<double> > transP(nStates, vector<double>(nStates));
  vector<vector<double> > updateTransP(nStates, vector<double>(nStates));
  double unif=log(1.0/nStates);
  for (size_t i=0;i<nStates;i++)
    {
      for (size_t j=0;j<nStates;j++)
        {
	  if(i==j)
            {
	      transP[i][j]= unif; //log(1 - std::exp(beta) );
            }
	  else
            {
	      if ( j==3)
                {
		  transP[i][j]= unif; // beta - log(nStates-1) ;
                }
	      else
                {
		  transP[i][j]= unif; // beta - log(nStates-1);                
                }
            }
        }
    }
  int maxCov=(int)mean/2*(maxState+1);
  vector<vector<double> > emisP(nStates, vector<double>(maxCov+1));
  vector<vector<double> > updateEmisP(nStates, vector<double>(maxCov+1));  

  //
  // Cap coverage where hmm does not bother calculating.
  //
  for (size_t c=0; c < covBins.size(); c++) {
    for (size_t b=0; b < covBins[c].size(); b++) {
      covBins[c][b] = min(covBins[c][b], maxCov);
    }
  }

  for (size_t i=0;i<nStates;i++){
    for (size_t j=0;j<maxCov;j++){
      if (model == POIS) {	
	emisP[i][j]=LgPrpoiss( (int) i , j , (int) mean/2 );
      }
      else  {
	emisP[i][j]=LgNegBinom((int)i, (int) j, mean/2, var/2);
      }      
    }
  }

  printModel(transP);
  printEmissions(emisP);
  
  //
  // This data will be used for all threads.
  //
  
  for (int procIndex = 0; procIndex < nproc; procIndex++) {
    pthread_attr_init(&threadAttr[procIndex]);
    if (bamFileName != "") {
      threadInfo[procIndex].htsfp = hts_open(bamFileName.c_str(),"r");
      
    }
    else {
      threadInfo[procIndex].htsfp = NULL;
    }
    threadInfo[procIndex].bamidx = bamidx;
    threadInfo[procIndex].samHeader=samHeader;
    threadInfo[procIndex].fai = fai_load_format(referenceName.c_str(), FAI_FASTA);
    
    threadInfo[procIndex].lastSeq = &curSeq;
    threadInfo[procIndex].semaphore = &semaphore;
    threadInfo[procIndex].contigNames = &contigNames;
    threadInfo[procIndex].contigLengths = &contigLengths;
    threadInfo[procIndex].covBins = &covBins;
    threadInfo[procIndex].copyNumber = &copyNumber;    
    threadInfo[procIndex].snvs = &snvs;    
    threadInfo[procIndex].lepsi = lepsi;
    threadInfo[procIndex].scale = scale;    
    threadInfo[procIndex].maxState = maxState;
    threadInfo[procIndex].maxCov = maxCov;
    threadInfo[procIndex].mean = mean;
    threadInfo[procIndex].var = var;
    threadInfo[procIndex].transP = &transP;
    threadInfo[procIndex].emisP = &emisP;
    threadInfo[procIndex].startP = &startP;
  }

  if (covBedInFileName == "") {
    for (int procIndex = 0; procIndex < nproc; procIndex++) {  
      pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ParseChrom, &threadInfo[procIndex]);
    }
      
    for (int procIndex = 0; procIndex < nproc; procIndex++) {
      pthread_join(threads[procIndex], NULL);
    }

    if (covBedOutFileName != "" ) {
      WriteCovBed(covBedOutFileName, contigNames, covBins);
    }
  }

  
  printModel(transP);
  for (int i=0; i < 10; i++) {
    BaumWelchEMStep(startP, transP, emisP,
		    covBins[0],
		    f[0],b[0],
		    updateTransP, updateEmisP);
    
    cout << "update" <<endl;
    printModel(updateTransP);
    transP=updateTransP;
    //
    // Gather per-state summary statistics
    //
    vector<long> stateTotCov(nStates, 0);
    vector<int>  stateNSamples(nStates, 0);
    
    //emisP=updateEmisP;
    
    printModel(transP);
    Reset(f);
    Reset(b);
  }
  /*
  ofstream outFile;
  ostream *outPtr;
  if (outFileName == "" or outFileName == "-" or outFileName == "stdin" or outFileName == "/dev/stdin") {
    outPtr = &cout;
  }
  else {
    outFile.open(outFileName.c_str());
      outPtr = &outFile;
  }
  */

  /*
  if (mergeBins == false) {
    for (int c=0; c < contigNames.size(); c++) {
      for (int b=0; b < copyNumber[c].size(); b++) {
	(*outPtr) << contigNames[c] << "\t" << b*BIN_LENGTH << "\t" << min((b+1)*BIN_LENGTH, contigLengths[c]) << "\t" << covBins[c][b] << "\t" << copyNumber[c][b] << endl;
      }
    }
  }
  else {
    for (int c=0; c < contigNames.size(); c++) {
      int start=0;
      int b=0;
      while (b < copyNumber[c].size()) {
	long totCov=0;	
	while (b < copyNumber[c].size() and copyNumber[c][b] == copyNumber[c][start]) { totCov+=covBins[c][b]; b++;}
	(*outPtr) << contigNames[c] << "\t" << start*BIN_LENGTH << "\t" << b *BIN_LENGTH << "\t" << totCov/(b-start) << "\t" << copyNumber[c][start] << endl;
	start=b;
      }
    }
  }
  */
  /*

  int state;
  for (int c=0; c < contigNames.size(); c++) {
    for (int b=0; b < covBins[c].size(); b++) {
      state=ceil(covBins[c][b]/mean);
      maxState=max(state,maxState);
    }
  }
  */
  /*
  printModel(transP);
  printEmissions(emisP);
  */


}
