all: samToBed viterbi


viterbi: viterbi.cpp
	g++ -O2 -I $(CONDA_PREFIX)/include $^ -o $@

samToBed: SamToBed.cpp
	g++ $^ -o $@
