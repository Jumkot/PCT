#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N 20

struct particle {
    float x, y, z;
};

const float G = 6.67e-11;

////////////////////////// ПОСЛЕДОВАТЕЛЬНАЯ ВЕРСИЯ //////////////////////////

void calculate_forces(struct particle* p, struct particle* f, float* m)
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            // Вычисление силы, действующей на тело i со стороны j
            float dist = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2) + powf(p[i].z - p[j].z, 2));
            float mag = (G * m[i] * m[j]) / powf(dist, 2);
            struct particle dir = {
                    .x = p[j].x - p[i].x,
                    .y = p[j].y - p[i].y,
                    .z = p[j].z - p[i].z};
            // Сумма сил, действующих на тело i
            f[i].x += mag * dir.x / dist;
            f[i].y += mag * dir.y / dist;
            f[i].z += mag * dir.z / dist;
            // Сумма сил, действующих на тело j (симметричность)
            f[j].x -= mag * dir.x / dist;
            f[j].y -= mag * dir.y / dist;
            f[j].z -= mag * dir.z / dist;
        }
    }
}

void move_particles(struct particle* p, struct particle* f, struct particle* v, float* m, float dt)
{
    for (int i = 0; i < N; i++)
    {
        struct particle dv = {
                .x = f[i].x / m[i] * dt,
                .y = f[i].y / m[i] * dt,
                .z = f[i].z / m[i] * dt,
        };
        struct particle dp = {
                .x = (v[i].x + dv.x / 2) * dt,
                .y = (v[i].y + dv.y / 2) * dt,
                .z = (v[i].z + dv.z / 2) * dt,
        };
        v[i].x += dv.x;
        v[i].y += dv.y;
        v[i].z += dv.z;
        p[i].x += dp.x;
        p[i].y += dp.y;
        p[i].z += dp.z;
        f[i].x = f[i].y = f[i].z = 0;
    }
}

double nbody_serial()
{
    struct particle* p = malloc(sizeof(*p) * N); // Положение частиц (x, y, z)
    struct particle* f = malloc(sizeof(*f) * N); // Сила, действующая на каждую частицу (x, y, z)
    struct particle* v = malloc(sizeof(*v) * N); // Скорость частицы (x, y, z)
    float* m = malloc(sizeof(*m) * N); // Масса частицы
    for (int i = 0; i < N; i++)
    {
        p[i].x = rand() / (float)RAND_MAX - 0.5;
        p[i].y = rand() / (float)RAND_MAX - 0.5;
        p[i].z = rand() / (float)RAND_MAX - 0.5;
        v[i].x = rand() / (float)RAND_MAX - 0.5;
        v[i].y = rand() / (float)RAND_MAX - 0.5;
        v[i].z = rand() / (float)RAND_MAX - 0.5;
        m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
        f[i].x = f[i].y = f[i].z = 0;
    }
    double dt = 1e-5;

    double time = omp_get_wtime();
    
    for (double t = 0; t <= 1; t += dt)
    { // Цикл по времени (модельному)
        calculate_forces(p, f, m); // Вычисление сил – O(N^2)
        move_particles(p, f, v, m, dt); // Перемещение тел O(N)
    }

    time = omp_get_wtime() - time;

    free(m);
    free(v);
    free(f);
    free(p);

    return time;
}

////////////////////////// ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ //////////////////////////

void calculate_forces_parallel(struct particle* p, struct particle* f[], float* m)
{
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    for (int i = 0; i < N; i++) {
        f[tid][i].x = 0;
        f[tid][i].y = 0;
        f[tid][i].z = 0;
    }
    #pragma omp for schedule(dynamic, 8)
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            float dist = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2) + powf(p[i].z - p[j].z, 2));
            float mag = (G * m[i] * m[j]) / powf(dist, 2);
            struct particle dir = {
                    .x = p[j].x - p[i].x,
                    .y = p[j].y - p[i].y,
                    .z = p[j].z - p[i].z};
            f[tid][i].x += mag * dir.x / dist;
            f[tid][i].y += mag * dir.y / dist;
            f[tid][i].z += mag * dir.z / dist;
            f[tid][j].x -= mag * dir.x / dist;
            f[tid][j].y -= mag * dir.y / dist;
            f[tid][j].z -= mag * dir.z / dist;
        }
    }
    #pragma omp single // Итоговый вектор сил сформируем в первой строке – f[0][i]
    {
        for (int i = 0; i < N; i++)
        {
            for (int tid = 1; tid < nthreads; tid++)
            {
                f[0][i].x += f[tid][i].x;
                f[0][i].y += f[tid][i].y;
                f[0][i].z += f[tid][i].z;
            }
        }
    }
}

void move_particles_parallel(struct particle* p, struct particle* f[], struct particle* v, float* m, float dt)
{
    #pragma omp for
    for (int i = 0; i < N; i++)
    {
        struct particle dv = {
                .x = f[0][i].x / m[i] * dt,
                .y = f[0][i].y / m[i] * dt,
                .z = f[0][i].z / m[i] * dt,
        };
        struct particle dp = {
                .x = (v[i].x + dv.x / 2) * dt,
                .y = (v[i].y + dv.y / 2) * dt,
                .z = (v[i].z + dv.z / 2) * dt,
        };
        v[i].x += dv.x;
        v[i].y += dv.y;
        v[i].z += dv.z;
        p[i].x += dp.x;
        p[i].y += dp.y;
        p[i].z += dp.z;
    }
}

double nbody_parallel(int n_threads)
{
    struct particle* p = malloc(sizeof(*p) * N); // Положение частиц (x, y, z)
    struct particle* f[n_threads];
    for (int i = 0; i < n_threads; i++)
    {
        f[i] = malloc(sizeof(struct particle) * N);
        for (int j = 0; j < N; j++)
        {
            f[i][j].x = 0;
            f[i][j].y = 0;
            f[i][j].z = 0;
        }
    }
    struct particle* v = malloc(sizeof(*v) * N); // Скорость частицы (x, y, z)
    float* m = malloc(sizeof(*m) * N); // Масса частицы
    for (int i = 0; i < N; i++)
    {
        p[i].x = rand() / (float)RAND_MAX - 0.5;
        p[i].y = rand() / (float)RAND_MAX - 0.5;
        p[i].z = rand() / (float)RAND_MAX - 0.5;
        v[i].x = rand() / (float)RAND_MAX - 0.5;
        v[i].y = rand() / (float)RAND_MAX - 0.5;
        v[i].z = rand() / (float)RAND_MAX - 0.5;
        m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
    }
    double dt = 1e-5;

    double time = omp_get_wtime();
    
    #pragma omp parallel num_threads(n_threads)
    {
        for (double t = 0; t <= 1; t += dt)
        { // Цикл по времени (модельному)
            calculate_forces_parallel(p, f, m); // Вычисление сил – O(N^2)
            move_particles_parallel(p, f, v, m, dt); // Перемещение тел O(N)
        }
    }

    time = omp_get_wtime() - time;

    free(m);
    free(v);
    free(p);
    for (int i = 0; i < n_threads; i++)
    {
        free(f[i]);
    }

    return time;
}

int main()
{
    int t_serial = nbody_serial();

    FILE* file;
    file = fopen("data.dat", "w");
    for (int i = 2; i <= 8; i += 2)
    {
        fprintf(file, "%d    %f\n", i, t_serial / nbody_parallel(i));
    }
    fclose(file);

    return 0;
}