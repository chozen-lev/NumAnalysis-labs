/*
 *
 * Lab5, "СЕРЕДНЬОКВАДРАТИЧНЕ НАБЛИЖЕННЯ ФУНКЦІЙ"
 * Lev Bazylskyi, KV-51
 * 22. 10*x*x*cosh(x)*sin(13*x), [0;1.2]
 * (22 % 8 = 6 => 110). Многочлен Чебишева, схема з вибором головного елемента,
 * інтегрування коефіцієнтів нормальної системи за узагальненою формулою Сімпсона
*/

#include <iostream>
#include <math.h>
#include <fstream>
#include <cstring>

// Розкоментувати для виведення всіх результатів обчислень в консолі
// В файл, при цьому, нічого не записується
//#define DEBUG

#define A 0
#define B 1.2
#define EPS 1e-9 // максимальна похибка розрахунку інтеграла
#define N 26 // максимальна степінь полінома

#define DELTA_EPS 1e-2 // максимальне середньоквадратичне відхилення
#define NUM 100 // кількість точок для виведення в файл
#define FILE_PATH "table.csv"

#define f(x) (10*(x)*(x)*cosh(x)*sin(13*(x)))

double vector[N];

double phi(double x, int k)
{
    if (k == 0) {
        return 1.0;
    }
    if (k == 1) {
        return x;
    }

    double t0 = 1.0, t1 = x, buff;
    const double x2 = x + x;

    for (int i = 1; i < k; i++) {
        buff = x2*t1 - t0;
        t0 = t1;
        t1 = buff;
    }

    return t1;
}

double get_Pm(double x, int m, double vector[N])
{
    double result = 0.0;

    for (int i = 0; i < m; i++) {
        result += vector[i]*phi(x, i);
    }

    return result;
}

double Sigma(double x, int m, int)
{
    double result = f(x) - get_Pm(x, m, vector);

    return result*result;
}

double phi_product(double x, int i, int j)
{
    if (i == j) {
        double result = phi(x, i);
        return result*result;
    }
    if (j == N) {
        return f(x)*phi(x, i);
    }

    return phi(x, i)*phi(x, j);
}

double Simpson(double a, double b, double h, int n, double (*func)(double,int,int), int i, int j)
{
    double sigma1 = 0.0, sigma2 = 0.0;

    for (int k = 1; k < n; k += 2) {
        sigma1 += func(a + k*h, i, j);
    }

    for (int k = 2; k < n - 1; k += 2) {
        sigma2 += func(a + k*h, i, j);
    }

    return h/3*(func(a, i, j) + func(b, i, j) + 4*sigma1 + 2*sigma2);
}

double integral(double a, double b, double (*func)(double,int,int), int i, int j, double eps)
{
    long long int n = (b - a)/sqrt(sqrt(eps));
    double In = 0.0, I2n = 0.0, h = (b - a)/n;
    n = (b - a)/h;
    if(n % 2 == 1) {
        n++;
    }
    h = (b - a)/n;

    In = Simpson(a, b, h, n, func, i, j);
    n <<= 1; h /= 2;
    I2n = Simpson(a, b, h, n, func, i, j);

    eps *= 15;
    while (fabs((In - I2n)/I2n) > eps) {
        In = I2n;
        n <<= 1; h /= 2;
        I2n = Simpson(a, b, h, n, func, i, j);
    }
#ifdef DEBUG
    printf("[%d,%d] n = %lld, h = %f, value = %f\n", i, j, n, h, I2n);
#endif
    return I2n;
}

void get_matrix(double matrix[][N + 1], double a, double b, double eps)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            matrix[i][j] = integral(a, b, phi_product, i, j, eps);
        }
    }
}

void GaussianElimination(const double _matrix[][N + 1], double vector[N], const int n)
{
    double matrix[n][n + 1], max, buff;
    int maxRow;

    for (int i = 0; i < n; i++) {
        memcpy(matrix[i], _matrix[i], n*sizeof(double));
        matrix[i][n] = _matrix[i][N];
    }

    for (int k = 0; k < n; k++)
    {
        max = fabs(matrix[k][k]);
        maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(matrix[i][k]) > max) {
                max = fabs(matrix[i][k]);
                maxRow = i;
            }
        }

        for (int i = k; i < n + 1; i++) {
            buff = matrix[maxRow][i];
            matrix[maxRow][i] = matrix[k][i];
            matrix[k][i] = buff;
        }

        for (int i = k + 1; i < n; i++) {
            buff = -matrix[i][k]/matrix[k][k];
            for (int j = k; j < n + 1; j++) {
                if (k == j) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] += buff*matrix[k][j];
                }
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        vector[i] = matrix[i][n]/matrix[i][i];
        for (int k = i - 1; k >= 0; k--) {
            matrix[k][n] -= matrix[k][i]*vector[i];
        }
    }
}

double Delta(double a, double b, double eps, int m)
{
    return sqrt(integral(a, b, Sigma, m, -1, eps)/(b - a));
}

int main()
{
#ifdef DEBUG
    printf("Generalized polynomial:\n");
#else
    printf("Сalculation of a generalized polynomial ...\n");
#endif
    double matrix[N][N + 1];
    get_matrix(matrix, A, B, EPS);

#ifdef DEBUG
    printf("\nSearching polynomial degree where delta P(m) < %.0e\n", DELTA_EPS);
#else
    printf("Searching polynomial degree where delta P(m) < %.0e ...\n", DELTA_EPS);
#endif

    int i;
    double delta;
    for (i = 1; i < N; i++) {
        GaussianElimination(matrix, vector, i);
        delta = Delta(A, B, EPS, i + 1);
        if (delta < DELTA_EPS) {
        #ifdef DEBUG
            printf("m = %d\n", i + 1);
            printf("\nRoots of a generalized polynomial:\n");
            for (int j = 0; j < i + 1; j++) {
                printf("a%d = %e\n", j, vector[j]);
            }
        #endif
        #ifndef DEBUG
             std::ofstream file(FILE_PATH, std::ofstream::out);
             printf("Writing XY coordinates to file \"%s\" for delta P(%d) = %.2e ...\n\n", FILE_PATH, i + 1, delta);
        #else
             printf("\nXY coordinates:\n");
        #endif
            const double h = (double)(B - A)/NUM;
            for (double x = A; x <= B; x += h)
            {
            #ifdef DEBUG
                printf("%f, %f\n", x, get_Pm(x, N, vector));
            #else
                file << x << ";" << get_Pm(x, N, vector) << std::endl;
            #endif
            }
        #ifndef DEBUG
            file.close();
        #endif
            break;
        }
    }
    if (i == N) {
        printf("Something went wrong. Try to reduce DELTA_EPS or increase N.\n");
    }
    return 0;
}
