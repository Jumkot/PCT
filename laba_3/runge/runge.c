#include <stdio.h>
#include <math.h>
#include <omp.h>

#define EPS 1e-5
#define N 100000000

const double a = -1, b = 1;

double f(double x)
{
    return sin(x + 2) / (0.4 + cos(x));
}

void integral()
{
    int n = N, k;
    double sq[2], delta = 1;
    for (k = 0; delta > EPS; n *= 2, k ^= 1)
    {
        double h = (b - a) / n;
        double s = 0.0;
        for (int i = 0; i < n; i++)
        {
            s += f(a + h * (i + 0.5));
        }
        sq[k] = s * h;
        if (n > N)
        {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
}

void integral_omp(int n_threads)
{
    double sq[2];
    #pragma omp parallel num_threads(n_threads)
    {
        int n = N, k;
        double delta = 1;
        for (k = 0; delta > EPS; n *= 2, k ^= 1)
        {
            double h = (b - a) / n;
            double s = 0.0;
            sq[k] = 0;
            // Ждем пока все потоки закончат обнуление sq[k], s
            #pragma omp barrier

            #pragma omp for nowait
            for (int i = 0; i < n; i++)
            {
                s += f(a + h * (i + 0.5));
            }
            #pragma omp atomic
            sq[k] += s * h;
            // Ждем пока все потоки обновят sq[k]
            #pragma omp barrier

            if (n > N)
            {
                delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
            }
            #if 0
            printf("n=%ld i=%ld sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
            #endif
        }
    }
}

double run_serial()
{
    double t = omp_get_wtime();
    integral();
    t = omp_get_wtime() - t;
    printf("Elapsed time (serial): %.6f sec.\n", t);
    
    return t;
}

double run_parallel(int n_threads)
{
    double t = omp_get_wtime();
    integral_omp(n_threads);
    t = omp_get_wtime() - t;
    printf("Elapsed time (parallel: %d threads): %.6f sec.\n", n_threads, t);
    
    return t;
}

int main()
{
    double t_serial= run_serial();

    FILE* file;
    file = fopen("data.dat", "w");
    for (int i = 2; i <= 8; i += 2)
    {
        fprintf(file, "%d    %f\n", i, t_serial/ run_parallel(i));
    }
    fclose(file);

    return 0;
}