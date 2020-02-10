# hmm
HMM tools
 pipeline for calling CNVs in assmeblies or alignments
 
 needs :
 bedtools
 samtools
 python 3.6.7
 snakemake
 boost
 R
 
 
 viterbi.cpp needs boost to compile
 using g++ 7.3.0
  g++ -W -I {boost_install}/boost/include/ viterbi.cpp -o viterbi
  usage: ./viterbi <coverage observation> <Mean coverage> <output prefix> <scaler>
