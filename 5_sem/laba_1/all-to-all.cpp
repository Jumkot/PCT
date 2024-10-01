#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int message_size = 1024;

    std::vector<char> send_buffer(msg_size * size, rank);
    std::vector<char> recv_buffer(msg_size * size);

    std::vector<MPI_Request> requests(2 * size);

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    for (int i = 0; i < size; ++i) {
        MPI_Isend(
                &send_buffer[i * message_size],
                message_size,
                MPI_CHAR,
                i,
                0,
                MPI_COMM_WORLD,
                &requests[i]);
        MPI_Irecv(
                &recv_buffer[i * message_size],
                message_size,
                MPI_CHAR,
                i,
                0,
                MPI_COMM_WORLD,
                &requests[size + i]);
    }

    MPI_Waitall(2 * size, requests.data(), MPI_STATUSES_IGNORE);

    double end = MPI_Wtime();
    std::cout << "Rank: " << rank << ", Time: " << end - start << " seconds\n";

    MPI_Finalize();
    return 0;
}
