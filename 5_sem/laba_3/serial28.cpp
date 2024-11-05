#include <cmath>
#include <ctime>
#include <iostream>
#include <mpi.h>
#include <vector>

int n = 5000;

// Последовательная версия программы
void gauss_serial()
{
    std::vector<std::vector<double>> a(n, std::vector<double>(n + 1));
    std::vector<double> x(n, 0.0);

    // Инициализация матрицы и столбца свободных членов
    for (int i = 0; i < n; i++) {
        srand(i * (n + 1));
        for (int j = 0; j < n; j++) {
            a[i][j] = rand() % 100 + 1;
        }
        // Последний столбец (n-й) — свободные члены
        a[i][n] = rand() % 100 + 1;
    }

    // Прямой ход (исключение переменных) — сложность O(n^3)
    for (int k = 0; k < n - 1; k++) {
        double pivot = a[k][k];
        for (int i = k + 1; i < n; i++) {
            double lik = a[i][k] / pivot;
            for (int j = k; j < n + 1; j++) {
                a[i][j] -= lik * a[k][j];
            }
        }
    }

    // Обратный ход (вычисление переменных) — сложность O(n^2)
    for (int k = n - 1; k >= 0; k--) {
        x[k] = a[k][n];
        for (int i = k + 1; i < n; i++) {
            x[k] -= a[k][i] * x[i];
        }
        x[k] /= a[k][k];
    }
}

int main(int argc, char* argv[])
{
    std::cout.setf(std::ios::fixed);

    // Запуск последовательной версии и замер времени выполнения
    double serial_time = MPI_Wtime();
    gauss_serial();
    serial_time = MPI_Wtime() - serial_time;

    std::cout << serial_time << "\n";

    return 0;
}