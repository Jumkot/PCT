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
	sbatch job/*.job

sbatch_all:
	sbatch job/1000/*.job
	sbatch job/2500/*.job
	sbatch job/4000/*.job

out:
	g++ -o output output.cpp
	./output
	gnuplot graph.gp

.PHONY: all clean sed_all sbatch_all
