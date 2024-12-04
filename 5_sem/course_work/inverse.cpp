#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int rank, commsize, lb, ub, nrows;
int n;

// Вспомогательный класс для поиска максимума
class Data
{
public:
    double value;
    int index;

    Data(double value_, int index_): value(value_), index(index_) {}
    Data()
    {
        value = 0.0;
        index = 0;
    }    
};

// Вычисляет диапазон строк для текущего процесса
void get_chunk(int *u, int *l) {
    int rows = n / commsize; // Число строк на процесс
    int mod = n % commsize;  // Оставшиеся строки
    if (rank < mod) {
        *l = rank * (rows + 1);
        *u = *l + rows;
    } else {
        *l = mod * (rows + 1) + (rank - mod) * rows;
        *u = *l + rows - 1;
    }
    if (*u >= n) {
        *u = n - 1; // Гарантия, что верхний индекс не выходит за пределы
    }
}

// Определяет, какой процесс обрабатывает строку с данным индексом
int get_proc(int index) {
    int rows = n / commsize; // Число строк на процесс
    int mod = n % commsize;  // Оставшиеся строки
    int threshold = mod * (rows + 1); // Порог между процессами
    if (index < threshold) {
        return index / (rows + 1);
    } else {
        return mod + (index - threshold) / rows;
    }
}

// Создаёт матрицу заданного типа (диагональная или стартовая)
void create_matrix(double* matrix, bool base) {
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

// Выполняет шаг инверсии методом Гаусса-Жордана
// Все операции дублируются на матрице inverse
bool inverse_matrix(double* start, double* inverse, int current_col) {
    double local_max = 0.0;
    int local_index = -1;

    // Поиск максимального элемента в текущем столбце
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

    // Сравнение локальных максимумов для поиска глобального
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    if (global.value == 0.0) {
        return false; // Матрица необратима
    }

    int diagonal = get_proc(current_col); // Процесс с диагональным элементом
    int main = get_proc(global.index);    // Процесс с максимальным элементом

    double* local_buf = new double[4 * n];
    double* reduce_buf = new double[4 * n];
    std::fill(local_buf, local_buf + 4 * n, 0.0);
    std::fill(reduce_buf, reduce_buf + 4 * n, 0.0);

    // Сохранение строки с диагональным элементом
    if (rank == diagonal) {
        int index = current_col - lb; 
        for (int i = 0; i < n; i++) {
            local_buf[n + i] = start[index * n + i];
            local_buf[3 * n + i] = inverse[index * n + i];
        }
    }

    // Сохранение строки с максимальным элементом
    if (rank == main) {
        int index = global.index - lb;
        for (int i = 0; i < n; i++) {
            local_buf[i] = start[index * n + i];
            local_buf[2 * n + i] = inverse[index * n + i];
        }
    }

    // Сбор строк с макс. эл-тами и диаг. строк со всех процессов
    MPI_Allreduce(local_buf, reduce_buf, 4 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double c = reduce_buf[current_col]; // Значение диагонального элемента
    for (int i = 0; i < n; i++) {
        reduce_buf[i] /= c;             // Нормализация строки
        reduce_buf[2 * n + i] /= c;
    }

    // Обновление строки с максимальным элементом
    if (rank == main) {
        int index = global.index - lb;
        for (int i = 0; i < n; i++) {
            start[index * n + i] = reduce_buf[n + i];
            inverse[index * n + i] = reduce_buf[3 * n + i];
        }
    }

    // Обновление диагональной строки
    if (rank == diagonal) {
        int index = current_col - lb;
        for (int i = 0; i < n; i++) {
            start[index * n + i] = reduce_buf[i];
            inverse[index * n + i] = reduce_buf[2 * n + i];
        }
    }

    // Обновление остальных строк
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

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 1) {
        n = std::atoi(argv[1]);
    }

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

    // Сбор частей обратной матрицы на процессе 0
    double* full_inverse = nullptr;
    if (rank == 0) {
        full_inverse = new double[n * n]; // Полный массив для хранения данных
    }

    // Подготовка размеров и смещений для сбора матрицы на проццессе 0
    int* elem_per_process = new int[commsize];
    int* displace_per_process = new int[commsize];

    for (int i = 0; i < commsize; ++i) {
        int l, u;
        get_chunk(&u, &l);
        elem_per_process[i] = (u - l + 1) * n; // Количество элементов для каждого процесса
        displace_per_process[i] = l * n;              // Смещение в общем массиве
    }

    MPI_Gatherv(
        inverse, nrows * n, MPI_DOUBLE,        // Данные текущего процесса
        full_inverse, elem_per_process, displace_per_process,     // Данные всех процессов
        MPI_DOUBLE, 0, MPI_COMM_WORLD         // Тип и root
    );

    delete[] elem_per_process;
    delete[] displace_per_process;

    time += MPI_Wtime();

    if (rank == 0) {
        std::cout << commsize << "    " << time << "\n";
        delete[] full_inverse; // Освобождаем память, если она была выделена
    }

    delete[] start;
    delete[] inverse;

    MPI_Finalize();
    return 0;
}

