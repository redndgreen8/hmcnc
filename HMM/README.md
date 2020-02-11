# CNV hmm
## HMM tools
### pipeline for calling CNVs in assemblies or alignments
 
 needs :
 1. bedtools
 2. samtools
 3. python 3.6.7
 4. snakemake
 5. boost
 6. R
7. gcc 7.3.0

SamToBed can be compiled:
    `g++ -02 SamToBed.cpp -o samToBed`
 
viterbi.cpp needs boost to compile

    `g++ -W -I {boost_install}/boost/include/ viterbi.cpp -o viterbi`
    
viterbi coverage_observation Mean_coverage output_prefix scaler

bash hmm_caller.vert.sh bam_file ref.fai_file 0/1(for filtering clr reads)
