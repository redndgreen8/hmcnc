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
double lepsi=-500;

class SNV {
public:
  int ref;
  int alt;
  int pos;
};

using namespace std;
class ThreadInfo {
public:
  htsFile *htsfp;
  hts_idx_t *bamidx;
  bam_hdr_t *samHeader;			
  faidx_t *fai;
  int *lastSeq;
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


static void printModel(vector<vector<double> > &transP)
{

  cerr << "\nTRANS: \n";
  for (int r=0;r<transP.size();r++)
    {
      cerr << r <<": ";
      for (int c=0;c<transP[r].size();c++)
        {
	  cerr << transP[r][c] << " ";
        }
      cerr << "\n";
    }
  cerr<< "\n";
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



void viterbi( vector<double> &startP,
	      vector<vector<double> > &transP,
	      vector<vector<double> > &emisP,
	      vector<int> &observations,
	      size_t mean,
	      vector<int> &  viterbiPath, int maxAllowedCov){

  //size_t  nObservations  = observations.size();
  size_t nObservations=observations.size();
  int nStates=transP.size();
  vector<vector<double> > v(nStates, vector<double>(observations.size())  );
  vector<vector<double> > opt(nStates, vector<double>(observations.size())  );  
  // Init
  int obs = std::min(maxAllowedCov , observations[0]);
  for(size_t i=0;i<nStates;i++)
    {
      v[i][0] =  startP[i] + emisP[i][obs] ;
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
  /*
  for (size_t k=1; k <nObservations; k++) {
    cout << k << "\t" << observations[k] << "\t";
    for (size_t i=0; i < nStates; i++) {
      cout << std::setw(8) << v[i][k] << " ";
    }
    cout << endl;
  }

  for (size_t k=1; k <nObservations; k++) {
    cout << k << "\t" << observations[k] << "\t";
    for (size_t i=0; i < n States; i++) {
      cout << std::setw(4 ) << opt[i][k] << " ";
    }
    cout << endl;
  }
  */
      
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
  //			cout << bam_get_qname(b)  << "\t" << refLen << "\t" << (int) b->core.qual << "\t" << (int) b->core.flag << endl;			
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
    
    bool continueParsing=true;
    vector<bam1_t*> reads; //(bam_init1());
    long totalSize=0;
    int chunkNumber=0;
    int totalReads=0;
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
	reads.push_back(b);
	++totalReads;
      }
      cerr << "Reading contig index " << curSeq << " (chunk " << chunkNumber << ")\t" << (*threadInfo->contigNames)[curSeq] << " " << reads.size() << "\t" << totalReads << endl;
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
    
    if ((*threadInfo->covBins)[curSeq].size() == 0) { return ;}
    
    //
    double chromMean;
    long chromTot=0;
    if ((*threadInfo->covBins)[curSeq].size() == 0) {
      continue;
    }    
    viterbi( *threadInfo->startP, *threadInfo->transP, *threadInfo->emisP, (*threadInfo->covBins)[curSeq], threadInfo->mean, (*threadInfo->copyNumber)[curSeq], threadInfo->maxCov);
  }
  pthread_exit(NULL);    
}
    


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
int EstimateCoverage(string &bamFileName, vector<string> &chroms, vector<int> &lengths, string &useChrom, double &mean, double &var) {

  if (useChrom == "") {
    int maxLen=0;
    for (int i=0; i < lengths.size(); i++) {
      if (lengths[i] > maxLen) {
	useChrom = chroms[i];
      }
    }
  }
  int contigLength=0;
  for (int i=0; i < chroms.size(); i++) {
    if (chroms[i] == useChrom) {
      contigLength=lengths[i];
      break;
    }
  }
  if (contigLength == 0) {
    cerr << "ERROR Could not estimate coverage." <<endl;
    exit(1);
  }
  
  htsFile *htsfp;
    
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
  vector<int> covBins(contigLength/100+1);
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
      covBins[binIndex] = binTot/100;

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
    if (nSamples > 40000) {
      return 1;
    }
  }
  return 1;
}






int main(int argc, const char* argv[]) {
  int nproc=4;
  double scale=10;

  if (argc < 3) {
    cout << "usage: hmmcnc input.bam reference.fa" << endl	       
	 << "    The per nucleotide frequency will be calculated " << endl
	 << "    for the bam file for every region specified in regions.txt" << endl
	 << "    This file has one region per line, in the format chrom:start-end" << endl
	 << " -t value (int) Number of threads (4) " << endl
	 << " -c contig  Use this contig to estimate coverage. By default, longest contig." << endl
	 << " -e value (float) Value of log-epsilon (-500)." << endl      
	 << " -s value (float) Scalar for transition probabilities (pre Baum-Welch) (10)" << endl
	 << " -m value [pois|nb] Coverage model to use, Poisson (pois), or negative binomial (nb). Default nb." << endl
	 << " -x value Max state to allow (10)" << endl;
    exit(1);
  }
  NucMap.resize(256,4);
  NucMap[(int)'A']=0;
  NucMap[(int)'C']=1;
  NucMap[(int)'G']=2;
  NucMap[(int)'T']=3;	
  int maxState=10;	
  string bamFileName=argv[1];
  string referenceName=argv[2];
  typedef enum { POIS, NEG_BINOM  } MODEL_TYPE;
  MODEL_TYPE model=NEG_BINOM;

  if (argc > 3) {
    int argi=3;
    while (argi < argc) {
      if (strcmp(argv[argi], "-t") == 0) {
	++argi;
	nproc=atoi(argv[argi]);
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
      ++argi;
    }
  }




  
  string faiFileName=referenceName + ".fai";
  ifstream faiIn(faiFileName.c_str());
  if (faiIn.good() == false) {
    cerr << "ERROR. Reference is not indexed, or could not open .fai file" << endl;
    exit(1);
  }

  vector<string> contigNames;
  vector<int>    contigLengths;
  while (faiIn) {
    string line;
    getline(faiIn, line);
    stringstream strm(line);
    if (line != "") {
      string contig;
      strm >> contig;
      contigNames.push_back(contig);
      int length;
      strm >> length;
      contigLengths.push_back(length);
    }
  }

  string empty="";
  double mean;
  double var;
  EstimateCoverage(bamFileName, contigNames, contigLengths, empty, mean, var);
  cerr << "Got cov " << mean << "\t" << var << endl;
  
  //
  // Get the header of the bam file
  //
  htsFile *htsfp;

  htsfp = hts_open(bamFileName.c_str(),"r");
  bam_hdr_t *samHeader;			
  samHeader = sam_hdr_read(htsfp);


  faidx_t *fai = fai_load_format(referenceName.c_str(), FAI_FASTA);

  // 
  // Read index for random io
  //
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

  pthread_t *threads = new pthread_t[nproc];
  vector<ThreadInfo> threadInfo(nproc);		
  pthread_mutex_t semaphore;		
  pthread_mutex_init(&semaphore, NULL);
  pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
  int curSeq=0;
  vector<vector<int> > covBins;
  vector<vector<SNV> > snvs;
  vector<vector<int > > copyNumber;  
  copyNumber.resize(contigLengths.size());

  covBins.resize(contigLengths.size());
  for (int c=0; c < contigLengths.size(); c++ ) {
    covBins[c].resize(contigLengths[c]/BIN_LENGTH);
    copyNumber[c].resize(contigLengths[c]/BIN_LENGTH);       
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

  // trans prob, scale 3->2 by overlap of pdf.

  poisson distribution1(3*mean/2);
  double result3=pdf(distribution1, 3*mean/2);

  poisson distribution2(2*mean/2);
  double result2=pdf(distribution2, 3*mean/2);

  double epsi23 = result2/result3;

  //*300

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
  int maxCov=(int)mean/2*(maxState+1);
  vector<vector<double> > emisP(nStates, vector<double>(maxCov));


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

  
  for (int procIndex = 0; procIndex < nproc; procIndex++) {
    pthread_attr_init(&threadAttr[procIndex]);
    threadInfo[procIndex].htsfp = hts_open(bamFileName.c_str(),"r");
    threadInfo[procIndex].bamidx = bamidx;
    threadInfo[procIndex].samHeader=samHeader;
    threadInfo[procIndex].fai = fai;
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

    pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ParseChrom, &threadInfo[procIndex]);
  }

  
  for (int procIndex = 0; procIndex < nproc; procIndex++) {
    pthread_join(threads[procIndex], NULL);
  }
  for (int c=0; c < contigNames.size(); c++) {
    for (int b=0; b < copyNumber[c].size(); b++) {
      cout << contigNames[c] << b*BIN_LENGTH << "\t" << min((b+1)*BIN_LENGTH, contigLengths[c]) << "\t" << covBins[c][b] << "\t" << copyNumber[c][b] << endl;
    }
  }
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
