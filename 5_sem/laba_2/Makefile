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

sed_all:
	sed -i 's/\r//' midpoint_4.job
	sed -i 's/\r//' midpoint_16.job
	sed -i 's/\r//' midpoint_28.job
	sed -i 's/\r//' midpoint_40.job
	sed -i 's/\r//' midpoint_52.job
	sed -i 's/\r//' midpoint_64.job
	sed -i 's/\r//' monte-carlo_4.job
	sed -i 's/\r//' monte-carlo_16.job
	sed -i 's/\r//' monte-carlo_28.job
	sed -i 's/\r//' monte-carlo_40.job
	sed -i 's/\r//' monte-carlo_52.job
	sed -i 's/\r//' monte-carlo_64.job

sbatch_all:
	sbatch ./midpoint_4.job
	sbatch ./midpoint_16.job
	sbatch ./midpoint_28.job
	sbatch ./midpoint_40.job
	sbatch ./midpoint_52.job
	sbatch ./midpoint_64.job
	sbatch ./monte-carlo_4.job
	sbatch ./monte-carlo_16.job
	sbatch ./monte-carlo_28.job
	sbatch ./monte-carlo_40.job
	sbatch ./monte-carlo_52.job
	sbatch ./monte-carlo_64.job
	
cat_all:
	cat midpoint_4.txt
	cat midpoint_16.txt
	cat midpoint_28.txt
	cat midpoint_40.txt
	cat midpoint_52.txt
	cat midpoint_64.txt
	cat monte-carlo_4.txt
	cat monte-carlo_16.txt
	cat monte-carlo_28.txt
	cat monte-carlo_40.txt
	cat monte-carlo_52.txt
	cat monte-carlo_64.txt


.PHONY: all clean sed_all sbatch_all сat_all
