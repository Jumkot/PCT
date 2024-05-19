#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define THRESHOLD 1000
#define SIZE 100000

void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

void partition(int *v, int low, int high, int *left, int *right)
{
    int i = low;
    int j = high;
    int pivot = v[(i + j) / 2];

    while (i <= j)
    {
        while (v[i] < pivot)
        {
            i++;
        }
        while (v[j] > pivot)
        {
            j--;
        }

        if (i <= j)
        {
            swap(&v[i], &v[j]);
            i++;
            j--;
        }
    }

    *left = i;
    *right = j;
}

void quicksort_serial(int *v, int low, int high)
{
    int left = 0;
    int right = 0;
    partition(v, low, high, &left, &right);

    if (low < right)
    {
        quicksort_serial(v, low, right);
    }

    if (left < high)
    {
        quicksort_serial(v, left, high);
    }
}

void quicksort_tasks(int *v, int low, int high)
{
    int left = 0;
    int right = 0;
    partition(v, low, high, &left, &right);

    if (high - low < THRESHOLD || (right - low < THRESHOLD || high - left < THRESHOLD))
    {
        if (low < right)
        {
            quicksort_tasks(v, low, right);
        }
        if (left < high)
        {
            quicksort_tasks(v, left, high);
        }
    } else {
        #pragma omp task untied
        {
            quicksort_tasks(v, low, right);
        }
        quicksort_tasks(v, left, high);
    }
}

double run_parallel(int *array, int size, int n_threads)
{
    double t = omp_get_wtime();
    omp_set_num_threads(n_threads);

    #pragma omp parallel
    {
        #pragma omp single
        quicksort_tasks(array, 0, size - 1);
    }
    t = omp_get_wtime() - t;
    
    return t;
}

double run_serial(int *array, int size)
{
    double t = omp_get_wtime();
    quicksort_serial(array, 0, size - 1);    
    t = omp_get_wtime() - t;
    
    return t;
}

int main()
{
    char names[][10]   = { "data1.dat", "data2.dat", "data3.dat", "data4.dat", "data5.dat", "data6.dat" };

    for (int j = 1; j <= 5; j++)
    {
        int *array_serial = malloc(SIZE * j * sizeof(int));
        for (int i = 0; i < (SIZE * j); i++)
        {
            array_serial[i] = rand() % (SIZE * j) + 1;
        }

        FILE* file;
        char filename[10];
        strcpy(filename, names[j - 1]);
        file = fopen(filename, "w");

        double t_serial = run_serial(array_serial, SIZE * j);

        for (int i = 2; i <= 8; i += 2)
        {
            int *array_parallel = malloc(SIZE * j * sizeof(int));
            for (int k = 0; k < (SIZE * j); k++)
            {
                array_parallel[k] = rand() % (SIZE * j) + 1;
            }
            fprintf(file, "%d    %f\n", i,  t_serial / run_parallel(array_parallel, SIZE * j, i));
            free(array_parallel);
        }
        
        fclose(file);

        free(array_serial);
    }

    int *array_threshold = malloc(THRESHOLD * sizeof(int));
    for (int i = 0; i < THRESHOLD; i++)
    {
        array_threshold[i] = rand() % THRESHOLD + 1;
    }

    FILE* file;
    char filename[10];
    strcpy(filename, names[5]);
    file = fopen(filename, "w");

    double t_serial = run_serial(array_threshold, THRESHOLD);

    for (int i = 2; i <= 8; i += 2)
    {
        int *array_threshold_parallel = malloc(THRESHOLD * sizeof(int));
        for (int k = 0; k < THRESHOLD; k++)
        {
            array_threshold_parallel[k] = rand() % THRESHOLD + 1;
        }
        fprintf(file, "%d    %f\n", i,  t_serial / run_parallel(array_threshold_parallel, THRESHOLD, i));
        free(array_threshold_parallel);
    }

    fclose(file);

    free(array_threshold);

    return 0;
}