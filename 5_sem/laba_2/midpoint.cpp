#include <cmath>
#include <iostream>
#include <mpi.h>

double serial_time = 0.0;

const double eps = 1E-6;
const int n0 = 1000000;
const double a = 0.1;
const double b = 0.5;

double func(double x)
{
    return x / std::pow(std::sin(2 * x), 3);
}

// Последовательная версия программы
void integral_serial()
{
    int n = n0, k;
    double sq[2], delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1) {
        double h = (b - a) / n;
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += func(a + h * (i + 0.5));
        }

        sq[k] = s * h;
        if (n > n0) {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
}

int main(int argc, char** argv)
{
    std::cout.setf(std::ios::fixed);

    int commsize, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // Запуск последовательной версии и замер времени выполнения
        serial_time = MPI_Wtime();
        integral_serial();
        serial_time = MPI_Wtime() - serial_time;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int n = n0, k;
    double sq[2], delta = 1;

    double parallel_time = MPI_Wtime();
    for (k = 0; delta > eps; n *= 2, k ^= 1) {
        int points_per_proc = n / commsize;
        int lb = rank * points_per_proc;
        int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);
        double h = (b - a) / n;
        double s = 0.0;

        for (int i = lb; i <= ub; i++) {
            s += func(a + h * (i + 0.5));
        }

        MPI_Allreduce(&s, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sq[k] *= h;

        if (n > n0) {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
    parallel_time = MPI_Wtime() - parallel_time;

    if (rank == 0) {
        std::cout << commsize << "    " << serial_time / parallel_time << "\n";
    }

    MPI_Finalize();
    return 0;
}
