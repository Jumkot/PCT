CXX = mpicxx

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) -o $@ $<

clean:
	rm -f $(TARGETS)

sed_all:
	sed -i 's/\r//' job/gauss28_4.job
	sed -i 's/\r//' job/gauss28_8.job
	sed -i 's/\r//' job/gauss28_16.job
	sed -i 's/\r//' job/gauss28_32.job
	sed -i 's/\r//' job/gauss45_4.job
	sed -i 's/\r//' job/gauss45_8.job
	sed -i 's/\r//' job/gauss45_16.job
	sed -i 's/\r//' job/gauss45_32.job

sbatch_all:
	sbatch job/gauss28_4.job
	sbatch job/gauss28_8.job
	sbatch job/gauss28_16.job
	sbatch job/gauss28_32.job
	sbatch job/gauss45_4.job
	sbatch job/gauss45_8.job
	sbatch job/gauss45_16.job
	sbatch job/gauss45_32.job

out:
	g++ -o output output.cpp
	./output
	gnuplot graph.gp

.PHONY: all clean sed_all sbatch_all
