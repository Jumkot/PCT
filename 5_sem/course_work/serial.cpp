#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>

int n;

bool inverse_matrix(double* matrix) {
    double* common = new double[n * 2 * n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            common[i * 2 * n + j] = std::min(n - j, n - i);
        }
        for (int j = 0; j < n; ++j) {
            common[i * 2 * n + j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    for (int i = 0; i < n; ++i) {
        double diag = common[i * 2 * n + i];
        if (diag == 0.0) {
            std::cerr << "No inverse matrix\n";
            delete[] common;
            return false;
        }

        for (int j = 0; j < 2 * n; ++j) {
            common[i * 2 * n + j] /= diag;
        }

        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = common[k * 2 * n + i];
            for (int j = 0; j < 2 * n; ++j) {
                common[k * 2 * n + j] -= factor * common[i * 2 * n + j];
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i * n + j] = common[i * 2 * n + j + n];
        }
    }

    delete[] common;
    return true;
}

int main(int argc, char** argv) {
    std::cout.setf(std::ios::fixed);

    MPI_Init(&argc, &argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    n = std::atoi(argv[1]);

    double time = -MPI_Wtime();

    double* matrix = new double[n * n];

    if (!inverse_matrix(matrix)) {
        delete[] matrix;
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    time += MPI_Wtime();

    std::cout << time << "\n";

    delete[] matrix;
    MPI_Finalize();
    return 0;
}
