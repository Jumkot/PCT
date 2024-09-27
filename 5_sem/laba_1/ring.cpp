#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int message_size = 1; message_size <= 1024 * 1024;
         message_size *= 1024) {
        std::vector<char> message(message_size, rank);
        int prev = (rank - 1 + size) % size;
        int next = (rank + 1) % size;
        std::vector<char> recv_buffer(message_size);

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

        for (int i = 0; i < size - 1; ++i) {
            MPI_Sendrecv(
                    &message[0],
                    message_size,
                    MPI_CHAR,
                    next,
                    0,
                    &recv_buffer[0],
                    message_size,
                    MPI_CHAR,
                    prev,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);

            message = recv_buffer;
        }

        double end = MPI_Wtime();

        if (rank == 0) {
            std::cout << "Message size: " << message_size
                      << " bytes, Time: " << end - start << " seconds\n";
        }
    }

    MPI_Finalize();
    return 0;
}
