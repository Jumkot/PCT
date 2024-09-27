#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> msg_sizes = {1024, 1024 * 1024};

    for (int message_size = 1024; message_size <= 1024 * 1024;
         message_size *= 1024) {
        std::vector<char> buffer(message_size, rank);

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        MPI_Bcast(&buffer[0], message_size, MPI_CHAR, 0, MPI_COMM_WORLD);

        double end = MPI_Wtime();

        std::cout << "Rank: " << rank << ", Message size: " << message_size
                  << " bytes, Time: " << end - start << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
