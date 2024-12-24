#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mpi.h>

int rank, commsize, lb, ub, nrows;
int n;

class Data {
public:
    double value;
    int index;

    Data(double value_, int index_) : value(value_), index(index_)
    {
    }
    Data()
    {
        value = 0.0;
        index = 0;
    }
};

void get_chunk(int* u, int* l)
{
    int rows = n / commsize;
    int mod = n % commsize;
    if (rank < mod) {
        *l = rank * (rows + 1);
        *u = *l + rows;
    } else {
        *l = mod * (rows + 1) + (rank - mod) * rows;
        *u = *l + rows - 1;
    }
    if (*u >= n) {
        *u = n - 1;
    }
}

int get_proc(int index)
{
    int rows = n / commsize;
    int mod = n % commsize;
    int threshold = mod * (rows + 1);
    if (index < threshold) {
        return index / (rows + 1);
    } else {
        return mod + (index - threshold) / rows;
    }
}

void create_matrix(double* matrix, bool base)
{
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < n; ++j) {
            if (base) {
                matrix[i * n + j] = std::min(n - j, n - (i + lb));
            } else if (i + lb == j) {
                matrix[i * n + j] = 1.0;
            } else {
                matrix[i * n + j] = 0.0;
            }
        }
    }
}

bool inverse_matrix(double* start, double* inverse, int current_col)
{
    double local_max = 0.0;
    int local_index = -1;

    for (int i = 0; i < nrows; i++) {
        int global_index = i + lb;
        if (global_index >= current_col) {
            double value = std::fabs(start[i * n + current_col]);
            if (value > local_max) {
                local_max = value;
                local_index = global_index;
            }
        }
    }

    Data global(0.0, 0);
    Data local(local_max, local_index);

    MPI_Allreduce(
            &local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    if (global.value == 0.0) {
        return false;
    }

    int diagonal = get_proc(current_col);
    int main = get_proc(global.index); 

    double* local_buf = new double[4 * n];
    double* reduce_buf = new double[4 * n];
    std::fill(local_buf, local_buf + 4 * n, 0.0);
    std::fill(reduce_buf, reduce_buf + 4 * n, 0.0);

    if (rank == diagonal) {
        int index = current_col - lb;
        for (int i = 0; i < n; i++) {
            local_buf[n + i] = start[index * n + i];
            local_buf[3 * n + i] = inverse[index * n + i];
        }
    }

    if (rank == main) {
        int index = global.index - lb;
        for (int i = 0; i < n; i++) {
            local_buf[i] = start[index * n + i];
            local_buf[2 * n + i] = inverse[index * n + i];
        }
    }

    MPI_Allreduce(
            local_buf, reduce_buf, 4 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double c = reduce_buf[current_col];
    for (int i = 0; i < n; i++) {
        reduce_buf[i] /= c;
        reduce_buf[2 * n + i] /= c;
    }

    if (rank == main) {
        int index = global.index - lb;
        for (int i = 0; i < n; i++) {
            start[index * n + i] = reduce_buf[n + i];
            inverse[index * n + i] = reduce_buf[3 * n + i];
        }
    }

    if (rank == diagonal) {
        int index = current_col - lb;
        for (int i = 0; i < n; i++) {
            start[index * n + i] = reduce_buf[i];
            inverse[index * n + i] = reduce_buf[2 * n + i];
        }
    }

    for (int i = 0; i < nrows; i++) {
        if (i + lb != current_col) {
            c = start[i * n + current_col];
            for (int j = 0; j < n; j++) {
                start[i * n + j] -= reduce_buf[j] * c;
                inverse[i * n + j] -= reduce_buf[2 * n + j] * c;
            }
        }
    }

    delete[] local_buf;
    delete[] reduce_buf;

    return true;
}

int main(int argc, char** argv)
{
    std::cout.setf(std::ios::fixed);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n>\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    n = std::atoi(argv[1]);

    double time = -MPI_Wtime();

    get_chunk(&ub, &lb);
    nrows = ub - lb + 1;

    double* start = new double[nrows * n];
    double* inverse = new double[nrows * n];

    create_matrix(start, true);
    create_matrix(inverse, false);

    for (int i = 0; i < n; ++i) {
        if (!inverse_matrix(start, inverse, i)) {
            std::cerr << "No inverse matrix\n";
            delete[] start;
            delete[] inverse;
            MPI_Finalize();
            return 1;
        }
    }

    double* full_inverse = nullptr;
    if (rank == 0) {
        full_inverse = new double[n * n];
    }

    int* elem_per_process = new int[commsize];
    int* displace_per_process = new int[commsize];

    for (int i = 0; i < commsize; ++i) {
        int l, u;
        get_chunk(&u, &l);
        elem_per_process[i] = (u - l + 1) * n;
        displace_per_process[i] = l * n;
    }

    MPI_Gatherv(
            inverse,
            nrows * n,
            MPI_DOUBLE,
            full_inverse,
            elem_per_process,
            displace_per_process,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
    );

    time += MPI_Wtime();
    double global_time = 0;

    MPI_Reduce(&t, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        delete[] full_inverse;
    }
    
    std::cout << commsize << "    " << global_time << "\n";

    delete[] start;
    delete[] inverse;
    delete[] elem_per_process;
    delete[] displace_per_process;

    MPI_Finalize();
    return 0;
}
