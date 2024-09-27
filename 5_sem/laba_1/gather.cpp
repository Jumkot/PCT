#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int message_size = 1024; message_size <= 1024 * 1024;
         message_size *= 1024) {
        std::vector<char> send_buffer(message_size, rank);
        std::vector<char> recv_buffer;

        if (rank == 0) {
            recv_buffer.resize(message_size * size);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        MPI_Gather(
                &send_buffer[0],
                message_size,
                MPI_CHAR,
                recv_buffer.data(),
                message_size,
                MPI_CHAR,
                0,
                MPI_COMM_WORLD);

        double end = MPI_Wtime();

        if (rank == 0) {
            std::cout << "Message size: " << message_size
                      << " bytes, Time: " << end - start << " seconds\n";
        }
    }

    MPI_Finalize();
    return 0;
}
