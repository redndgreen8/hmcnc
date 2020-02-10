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


SamToBed can be compiled:
g++ -02 SamToBed.cpp -o samToBed
 
viterbi.cpp needs boost to compile

using g++ 7.3.0
`g++ -W -I {boost_install}/boost/include/ viterbi.cpp -o viterbi`
viterbi coverage_observation Mean_coverage output_prefix scaler
