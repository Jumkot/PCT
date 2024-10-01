CXX = mpicxx
CXXFLAGS = -Wall -Wextra -Werror

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)
	rm -f *.txt

run_ring_job_1:
	sed -i 's/\r//' ring_1.job
	sbatch ./ring_1.job

run_ring_job_2:
	sed -i 's/\r//' ring_2.job
	sbatch ./ring_2.job

run_broadcast_job:
	sed -i 's/\r//' broadcast.job
	sbatch ./broadcast.job

run_gather_job:
	sed -i 's/\r//' gather.job
	sbatch ./gather.job

run_all-to-all_job:
	sed -i 's/\r//' all-to-all.job
	sbatch ./all-to-all.job
	
jobs_cat:
	cat ring_1_file.txt
	cat ring_2_file.txt
	cat broadcast_file.txt
	cat gather_file.txt
	cat all-to-all_file.txt


.PHONY: all clean run_ring_job run_broadcast_job run_gather_job run_all-to-all_job
