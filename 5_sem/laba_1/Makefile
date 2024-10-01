CXX = mpicxx
CXXFLAGS = -Wall -Wextra -Werror

SOURCES = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCES))

all: $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)


#sbatch ./J_ring.job -> 882 (получение номера)
#squeue (проверка статуса)
#cat ./J_ring.o882 (вывод результата по имени задачи и номеру)

run_ring_job:
	sbatch ./ring.job

run_broadcast_job:
	sbatch ./broadcast.job

run_gather_job:
	sbatch ./gather.job

run_all-to-all_job:
	sbatch ./all-to-all.job

.PHONY: all clean run_ring_job run_broadcast_job run_gather_job run_all-to-all_job
