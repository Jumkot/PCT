CXX = mpicxx

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) -o $@ $<

clean:
	rm -f $(TARGETS)

sed_all:
	sed -i 's/\r//' job/laplace3_8.job
	sed -i 's/\r//' job/laplace3_16.job
	sed -i 's/\r//' job/laplace3_32.job
	sed -i 's/\r//' job/laplace4_8.job
	sed -i 's/\r//' job/laplace4_16.job
	sed -i 's/\r//' job/laplace4_32.job

sbatch_all:
	sbatch job/laplace3_8.job
	sbatch job/laplace3_16.job
	sbatch job/laplace3_32.job
	sbatch job/laplace4_8.job
	sbatch job/laplace4_16.job
	sbatch job/laplace4_32.job

out:
	g++ -o output output.cpp
	./output
	gnuplot graph.gp

.PHONY: all clean sed_all sbatch_all
