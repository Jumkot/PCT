#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

double getrand(unsigned int* seed)
{
    return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
    return x / (y * y);
}

void integral(int n)
{
    int in = 0;
    unsigned int seed = 0;
    double s = 0;
    for (int i = 0; i < n; i++) {
        double x = getrand(&seed);      /* x in [0, 1] */
        double y = getrand(&seed) * 5;  /* y in [2, 5] */
        if (y > 2)
        {
            in++;
            s += func(x, y);
        }
    }
    double v = in / n;
    double res = v * s / in;
    #if 0
    printf("Result: %.12f, n %d\n", res, n);
    #endif
}

void integral_omp(int n, int n_threads)
{
    int in = 0;
    double s = 0;
    #pragma omp parallel num_threads(n_threads)
    {
        double s_loc = 0;
        int in_loc = 0;
        unsigned int seed = omp_get_thread_num();
        #pragma omp for nowait
        for (int i = 0; i < n; i++) {
            double x = getrand(&seed);      /* x in [0, 1] */
            double y = getrand(&seed) * 5;  /* y in [2, 5] */
            if (y > 2)
            {
                in_loc++;
                s_loc += func(x, y);
            }
        }
        #pragma omp atomic
        s += s_loc;
        #pragma omp atomic
        in += in_loc;
    }
    double v = in / n;
    double res = v * s / in;
}

double run_serial(int n)
{
    double t = omp_get_wtime();
    integral(n);
    t = omp_get_wtime() - t;
    printf("Elapsed time (serial): %.6f sec.\n", t);
    
    return t;
}

double run_parallel(int n, int n_threads)
{
    double t = omp_get_wtime();
    integral_omp(n, n_threads);
    t = omp_get_wtime() - t;
    printf("Elapsed time (parallel: %d threads): %.6f sec.\n", n_threads, t);
    
    return t;
}

int main()
{
    int n = 10000000;
    printf("n = %d\n", n);
    double t_serial= run_serial(n);

    FILE* file7;
    file7 = fopen("data7.dat", "w");
    for (int i = 2; i <= 8; i += 2)
    {
        fprintf(file7, "%d    %f\n", i, t_serial/ run_parallel(n, i));
    }
    fclose(file7);

    n *= 10;
    printf("\nn = %d\n", n);
    t_serial= run_serial(n);

    FILE* file8;
    file8 = fopen("data8.dat", "w");
    for (int i = 2; i <= 8; i += 2)
    {
        fprintf(file8, "%d    %f\n", i, t_serial/ run_parallel(n, i));
    }
    fclose(file8);

    return 0;
}