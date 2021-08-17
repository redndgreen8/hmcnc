all: samToBed viterbi hmmcnc

hmmcnc: hmmcnc.cpp
	g++ -g -I $(CONDA_PREFIX)/include $^ -o $@ -L $(CONDA_PREFIX)/lib -lhts -lpthread

viterbi: viterbi.cpp
	g++ -O2 -I $(CONDA_PREFIX)/include $^ -o $@

samToBed: SamToBed.cpp
	g++ $^ -o $@
