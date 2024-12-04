CXX = mpicxx

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) -o $@ $<

clean:
	rm -f $(TARGETS)

sed_all:
	sed -i 's/\r//' job/*.job
	sed -i 's/\r//' job/1000/*.job
	sed -i 's/\r//' job/2500/*.job
	sed -i 's/\r//' job/4000/*.job

sbatch_serial:
	sbatch job/serial/serial1.job
	sbatch job/serial/serial2.job
	sbatch job/serial/serial4.job

sbatch_all:
	sbatch job/1000/inverse1_4.job
	sbatch job/1000/inverse1_8.job
	sbatch job/1000/inverse1_12.job
	sbatch job/1000/inverse1_16.job
	sbatch job/1000/inverse1_20.job
	sbatch job/1000/inverse1_24.job
	sbatch job/1000/inverse1_28.job
	sbatch job/1000/inverse1_32.job
	sbatch job/2500/inverse2_4.job
	sbatch job/2500/inverse2_8.job
	sbatch job/2500/inverse2_12.job
	sbatch job/2500/inverse2_16.job
	sbatch job/2500/inverse2_20.job
	sbatch job/2500/inverse2_24.job
	sbatch job/2500/inverse2_28.job
	sbatch job/2500/inverse2_32.job
	sbatch job/4000/inverse4_4.job
	sbatch job/4000/inverse4_8.job
	sbatch job/4000/inverse4_12.job
	sbatch job/4000/inverse4_16.job
	sbatch job/4000/inverse4_20.job
	sbatch job/4000/inverse4_24.job
	sbatch job/4000/inverse4_28.job
	sbatch job/4000/inverse4_32.job

out:
	g++ -o output output.cpp
	./output
	gnuplot graph.gp

.PHONY: all clean sed_all sbatch_all
