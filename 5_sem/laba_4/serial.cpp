#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>

#define EPS 0.001
#define PI 3.14159265358979323846
#define IND(i, j) ((i) * nx + (j))

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int n = std::atoi(argv[1]);
    int ny = n, nx = n;

    std::vector<double> local_grid(ny * nx, 0.0);
    std::vector<double> local_newgrid(ny * nx, 0.0);

    // Initialize top border: u(x, 0) = sin(pi * x)
    double dx = 1.0 / (nx - 1.0);
    for (int j = 0; j < nx; j++) {
        int ind = IND(0, j);
        local_newgrid[ind] = local_grid[ind] = sin(PI * dx * j);
    }
    // Initialize bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
    for (int j = 0; j < nx; j++) {
        int ind = IND(ny - 1, j);
        local_newgrid[ind] = local_grid[ind] = sin(PI * dx * j) * exp(-PI);
    }

    double ttotal = MPI_Wtime();

    int niters = 0;
    for (;;) {
        niters++;

        for (int i = 1; i < ny - 1; i++) { // Update interior points
            for (int j = 1; j < nx - 1; j++) {
                local_newgrid[IND(i, j)] =
                    (local_grid[IND(i - 1, j)] + local_grid[IND(i + 1, j)] +
                     local_grid[IND(i, j - 1)] + local_grid[IND(i, j + 1)]) * 0.25;
            }
        }
        // Check termination condition
        double maxdiff = 0;
        for (int i = 1; i < ny - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                int ind = IND(i, j);
                maxdiff = fmax(maxdiff, fabs(local_grid[ind] - local_newgrid[ind]));
            }
        }
        // Swap grids (after termination local_grid will contain result)
        std::swap(local_grid, local_newgrid);

        if (maxdiff < EPS)
            break;
    }

    ttotal = MPI_Wtime() - ttotal;

    std::cout << ttotal << "\n";

    MPI_Finalize();
    return 0;
}
