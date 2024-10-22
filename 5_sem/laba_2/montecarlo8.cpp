#include <cmath>
#include <ctime>
#include <iostream>
#include <mpi.h>

double serial_time = 0.0;

double get_rand()
{
    return static_cast<double>(rand()) / RAND_MAX;
}

double func(double x, double y)
{
    return std::pow(std::exp(x + y), 2);
}

// Последовательная версия программы
void integral_serial(int n)
{
    int in = 0;
    double s = 0;

    for (int i = 0; i < n; i++) {
        double x = get_rand();           // x in [0, 1]
        double y = get_rand() * (1 - x); // y in [0, 1 - x]

        in++;
        s += func(x, y);
    }

    double v = in / n;
    double res = v * s / in;
}

int main(int argc, char** argv)
{
    std::cout.setf(std::ios::fixed);

    int rank, commsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int n = 100000000;

    if (rank == 0) {
        // Запуск последовательной версии и замер времени выполнения
        serial_time = MPI_Wtime();
        integral_serial(n);
        serial_time = MPI_Wtime() - serial_time;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    srand(rank); // Уникальный seed для каждого процесса
    int in = 0;
    double s = 0;

    double parallel_start_time = MPI_Wtime();
    for (int i = rank; i < n; i += commsize) {
        double x = get_rand();           // x in [0, 1]
        double y = get_rand() * (1 - x); // y in [0, 1 - x]

        double f_value = func(x, y);

        in++;
        s += f_value;
    }

    // Сбор результатов от всех процессов
    int global_in = 0;
    MPI_Reduce(&in, &global_in, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double global_sum = 0.0;
    MPI_Reduce(&s, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double parallel_time = MPI_Wtime() - parallel_start_time;

    if (rank == 0) {
        double v = static_cast<double>(global_in) / n;
        double res = v * global_sum / global_in;
        std::cout << commsize << "    " << serial_time / parallel_time << "\n";
    }

    MPI_Finalize();
    return 0;
}
