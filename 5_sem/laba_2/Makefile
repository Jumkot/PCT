CXX = mpicxx

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) -o $@ $<

clean:
	rm -f $(TARGETS)
	rm -f *.txt
	rm -f *.dat

sed_all:
	sed -i 's/\r//' job/midpoint_4.job
	sed -i 's/\r//' job/midpoint_8.job
	sed -i 's/\r//' job/midpoint_16.job
	sed -i 's/\r//' job/midpoint_20.job
	sed -i 's/\r//' job/midpoint_28.job
	sed -i 's/\r//' job/midpoint_32.job
	sed -i 's/\r//' job/montecarlo7_4.job
	sed -i 's/\r//' job/montecarlo7_8.job
	sed -i 's/\r//' job/montecarlo7_16.job
	sed -i 's/\r//' job/montecarlo7_20.job
	sed -i 's/\r//' job/montecarlo7_28.job
	sed -i 's/\r//' job/montecarlo7_32.job
	sed -i 's/\r//' job/montecarlo8_4.job
	sed -i 's/\r//' job/montecarlo8_8.job
	sed -i 's/\r//' job/montecarlo8_16.job
	sed -i 's/\r//' job/montecarlo8_20.job
	sed -i 's/\r//' job/montecarlo8_28.job
	sed -i 's/\r//' job/montecarlo8_32.job

sbatch_all:
	sbatch job/midpoint_4.job
	sbatch job/midpoint_8.job
	sbatch job/midpoint_16.job
	sbatch job/midpoint_20.job
	sbatch job/midpoint_28.job
	sbatch job/midpoint_32.job
	sbatch job/montecarlo7_4.job
	sbatch job/montecarlo7_8.job
	sbatch job/montecarlo7_16.job
	sbatch job/montecarlo7_20.job
	sbatch job/montecarlo7_28.job
	sbatch job/montecarlo7_32.job
	sbatch job/montecarlo8_4.job
	sbatch job/montecarlo8_8.job
	sbatch job/montecarlo8_16.job
	sbatch job/montecarlo8_20.job
	sbatch job/montecarlo8_28.job
	sbatch job/montecarlo8_32.job

out:
	g++ -o output output.cpp
	./output
	gnuplot graph.gp

.PHONY: all clean sed_all sbatch_all
