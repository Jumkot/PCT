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


#sbatch ./J_ring.job -> 882 (получение номера)
#squeue (проверка статуса)
#cat ./J_ring.o882 (вывод результата по имени задачи и номеру)

run_ring_job:
	sed -i 's/\r//' ring.job
	sbatch ./ring.job
	cat ring_file.txt

run_broadcast_job:
	sed -i 's/\r//' broadcast.job
	sbatch ./broadcast.job
	cat gather_file.txt

run_gather_job:
	sed -i 's/\r//' gather.job
	sbatch ./gather.job
	cat ./broadcast

run_all-to-all_job:
	sed -i 's/\r//' all-to-all.job
	sbatch ./all-to-all.job
	cat all-to-all_file.txt

.PHONY: all clean run_ring_job run_broadcast_job run_gather_job run_all-to-all_job
