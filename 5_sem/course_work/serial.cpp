#include <iostream>
#include <vector>
#include <mpi.h>
#include <cmath>
#include <algorithm>

// Определяем тип Matrix как двумерный вектор
typedef std::vector<std::vector<double>> Matrix;

// Функция для инверсии матрицы
bool inverse_matrix(Matrix &matrix) {
    int n = matrix.size();
    Matrix temp(n, std::vector<double>(2 * n, 0.0));

    // Инициализация расширенной матрицы
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[i][j] = matrix[i][j];
        }
        temp[i][i + n] = 1.0;
    }

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        double diag = temp[i][i];
        if (diag == 0.0) // Матрица необратима
        {
            std::cerr << "No inverse matrix\n";
        }

        // Нормализация строки
        for (int j = 0; j < 2 * n; ++j) {
            temp[i][j] /= diag;
        }

        // Обнуление столбца
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = temp[k][i];
            for (int j = 0; j < 2 * n; ++j) {
                temp[k][j] -= factor * temp[i][j];
            }
        }
    }

    // Извлечение результата
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = temp[i][j + n];
        }
    }
    return true;
}

// Функция для создания матрицы
Matrix create_matrix(int n) {
    Matrix matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = std::min(n - j, n - i);
        }
    }
    return matrix;
}

int main(int argc, char **argv) {
    std::cout.setf(std::ios::fixed);

    MPI_Init(&argc, &argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int n = std::atoi(argv[1]);

    // Создание матрицы
    Matrix matrix = create_matrix(n);

    double time = MPI_Wtime();
    inverse_matrix(matrix);
    
    time = MPI_Wtime() - time;
    std::cout << time << "\n";

    MPI_Finalize();
    return 0;
}