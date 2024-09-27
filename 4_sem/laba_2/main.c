#include <stdio.h>
#include <omp.h>
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>

#define M 15000
#define N 15000
#define STEP 5000

void matrix_vector_product(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++) {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            c[i] += a[i * n + j] * b[j];
        }
    }
}

double run_serial(int m, int n)
{
    double *a, *b, *c;
    a = malloc(sizeof(*a) * m * n);
    b = malloc(sizeof(*b) * n);
    c = malloc(sizeof(*c) * m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = i + j;
        }
    }
    for (int j = 0; j < n; j++)
    {
        b[j] = j;
    }
    double t = omp_get_wtime();
    matrix_vector_product(a, b, c, m, n);
    t = omp_get_wtime() - t;
    printf("Elapsed time (serial): %.6f sec.\n\n", t);

    free(a);
    free(b);
    free(c);

    return t;
}

void matrix_vector_product_omp(double *a, double *b, double *c, int m, int n)
{
    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++) {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                c[i] += a[i * n + j] * b[j];
            }
        }
    }
}

double run_parallel(int m, int n, int num_of_threads)
{
    double *a = malloc(sizeof(*a) * m * n);
    double *b = malloc(sizeof(*b) * n);
    double *c = malloc(sizeof(*c) * m);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = i + j;
        }
    }
    for (int j = 0; j < n; j++)
    {
        b[j] = j;
    }

    double t = omp_get_wtime();
    omp_set_num_threads(num_of_threads);
    matrix_vector_product_omp(a, b, c, m, n);
    t = omp_get_wtime() - t;
    
    printf("Elapsed time (parallel): %.6f sec.\n", t);

    free(a);
    free(b);
    free(c);

    return t;
}

int main(int argc, char **argv)
{
    int m = M;
    int n = N;
    for (int j = 0; j < 3; j++)
    {
        printf("m = %d, n = %d\n", m, n);
        printf("Memory used: %" PRIu64 " MiB\n", ((m * n + m + n) * sizeof(double)) >> 20);
        double t = run_serial(m, n);
        for (int num_of_threads = 2; num_of_threads <= 8; num_of_threads++)
        {
            printf("Num of threads = %d\n", num_of_threads);
            double tp = run_parallel(m, n, num_of_threads);
            printf("Speedup = %.4lf\n\n", t / tp);
            num_of_threads++;
        }
        m += STEP;
        n += STEP;
    }

    return 0;
}