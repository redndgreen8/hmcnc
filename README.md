# markov CNV hunter
## HMM tools
### pipeline for calling CNVs in assemblies or alignments
 
Initially the required packages need to be installed.
On our linux cluster the easiest package management software is Anaconda/Miniconda. 

1. download shell script (64bit):
https://docs.conda.io/en/latest/miniconda.html#linux-installers

2. run script
`bash Miniconda3-latest-Linux-x86_64.sh`
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html


3. Project env - There are many ways to do this but you can set up a project specific environment with all the packages you need.
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands

bedtools
samtools
snakemake
boost
R 
gcc

`conda create --name <proj_env> bedtools samtools snakemake boost gcc=7.3.0`

conda install can be used to further add packages to <proj> environment with explicit version numbers.
 
`conda install -n <proj_env>  scipy=0.15.0`

always activate the env before attempting a run
`conda activate <proj_env>`

you might run into a conda init error the first time so run conda init and rerun


5.
SamToBed can be compiled: g++ -02 SamToBed.cpp -o samToBed

viterbi.cpp needs boost to compile

`g++ -W -I {boost_install}/boost/include/ viterbi.cpp -o viterbi`

6.
repeatMask.merged.bed has to be generated from repeat mask annotation in UCSC table browser then merge,
mergeBed -i <repeatmask.bed> > HMM/repeatMask.merged.bed

7. fai file has to be in the same dir as asm.fa

Finally can run the snakefile
`snakemake -p -s <script_dir(HMM)>/markovCNVhunter.snakefile -j <threads> --config subread=<0/1> asm=<hg38.fa/asm.fa> bam=<Genome.bam> t=<threads> MQ=<min_mapq_for_reads> script_dir=<script_dir(HMM)>`
