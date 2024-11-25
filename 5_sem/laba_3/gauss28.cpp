#include <cmath>
#include <ctime>
#include <iostream>
#include <mpi.h>
#include <vector>

double serial_time = 0.0;
int n = 28000;

int get_chunk(int commsize, int rank)
{
    int q = n / commsize;
    if (n % commsize) {
        q++;
    }
    int r = commsize * q - n;

    /* Compute chunk size for the process */
    int chunk = q;
    if (rank >= commsize - r)
    {
       chunk = q - 1;
    }
    return chunk;
}

int main(int argc, char* argv[])
{
    std::cout.setf(std::ios::fixed);
    
    int rank, commsize;
    MPI_Init(&argc, &argv);

    double parallel_start_time = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    
    int nrows = get_chunk(commsize, rank);
    
    std::vector<int> rows(nrows);  // Номера локальных строк
    std::vector<std::vector<double>> a(nrows, std::vector<double>(n + 1));
    // Матрица дополнена столбцом для вектора b
    std::vector<double> x(n, 0.0);
    std::vector<double> tmp(n + 1);

    // Инициализация как в последовательной версии
    for (int i = 0; i < nrows; i++) {
        rows[i] = rank + commsize * i;
        srand(rows[i] * (n + 1));
        for (int j = 0; j < n; j++)
        {
            a[i][j] = rand() % 100 + 1;
        }
        // b[i]
        a[i][n] = rand() % 100 + 1;
    }

    // Прямой ход
    int row = 0;
    for (int i = 0; i < n - 1; i++) {
        // Исключаем x_i
        if (i == rows[row]) {
            // Рассылаем строку i, находящуюся в памяти текущего процесса
            MPI_Bcast(a[row].data(), n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            tmp = a[row];
            row++;
        } else {
            MPI_Bcast(tmp.data(), n + 1, MPI_DOUBLE, i % commsize, MPI_COMM_WORLD);
        }
        // Вычитаем принятую строку из уравнений, хранящихся в текущем процессе
        if (row < nrows) {
            for (int j = row; j < nrows; j++) {
                double scaling = a[j][i] / tmp[i];
                for (int k = i; k < n + 1; k++)
                {
                    a[j][k] -= scaling * tmp[k];
                }
            }
        }
    }

    row = 0;
    for (int i = 0; i < n; i++) {
        if (row < nrows && i == rows[row]) {
            x[i] = a[row][n];
            row++;
        }
    }

    /* Обратный ход */
    row = nrows - 1;
    for (int i = n - 1; i > 0; i--) {
        if (row >= 0) {
            if (i == rows[row]) {
                x[i] /= a[row][i]; // Передаем найденное x_i
                MPI_Bcast(&x[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
                row--;
            } else {
                MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % commsize, MPI_COMM_WORLD);
            }
        } else {
            MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % commsize, MPI_COMM_WORLD);
        }
        for (int j = 0; j <= row; j++) // Корректировка локальных x_i
        {
            x[rows[j]] -= a[j][i] * x[i];
        }
    }
    if (rank == 0)
    {
        x[0] /= a[row][0]; // Корректировка x_0
    }
    MPI_Bcast(&x[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Все процессы содержат корректный вектор x[n] решений

    double parallel_time = MPI_Wtime() - parallel_start_time;
    if (rank == 0) {
        std::cout << commsize << "    " << parallel_time << "\n";
    }

    MPI_Finalize();

    return 0;
}