/*
 *
 * Lab6, "НАБЛИЖЕННЯ ФУНКЦІЙ ЗА ДОПОМОГОЮ ІНТЕРПОЛЯЦІЙНИХ СПЛАЙНІВ"
 * Lev Bazylskyi, KV-51
 * 22. 10*x*x*cosh(x)*sin(13*x), [0;1.2]
 * (22 % 6 = 4 => 100). Метод прогону
*/

#include <iostream>
#include <math.h>
#include <cstring>
#include <fstream>

#define A 0
#define B 1.2
#define N 100

#define FILE_PATH "table.csv"

#define f(x) (10*x*x*cosh(x)*sin(13*x))

typedef struct
{
    double a, b, c, d, x;
} Spline;

void get_matrix(double matrix[][N + 1], double a, double b)
{
    const double h = (b - a)/N;
    double x = a;

    memset(matrix, 0, N*(N + 1)*sizeof(double));

    for (int i = 1; i < N - 1; i++, x += h) {
        matrix[i][i - 1] = h; // A
        matrix[i][i] = 4*h; // C
        matrix[i][i + 1] = h; // B
        matrix[i][N] = 6*(f(x - h) - 2*f(x) + f(x + h))/h; // F
    }
    matrix[0][0] = 1.0;
    matrix[N - 1][N - 1] = 1.0;
}

// метод прогону
void Tridiagonal(double matrix[][N + 1], double vector[N])
{
    double alfa[N], beta[N];

    alfa[0] = -matrix[0][1]/matrix[0][0];
    beta[0] = matrix[0][N]/matrix[0][0];
    for (int i = 0; i < N - 1; i++) {
        alfa[i + 1] = -matrix[i][i + 1]/(matrix[i][i - 1]*alfa[i] + matrix[i][i]);
        beta[i + 1] = (matrix[i][N] - matrix[i][i - 1]*beta[i])/(matrix[i][i - 1]*alfa[i] + matrix[i][i]);
    }

    vector[N - 1] = (matrix[N - 1][N] - matrix[N - 2][N - 3]*beta[N - 1])/(matrix[N - 1][N - 2]*alfa[N - 1] + matrix[N - 1][N - 1]);
    for (int i = N - 2; i >= 0; i--) {
        vector[i] = alfa[i + 1]*vector[i + 1] + beta[i + 1];
    }
}

void get_splines(double a, double b, const double c[N], Spline splines[N])
{
    const double h = (b - a)/N;
    double x = a;

    splines[0].a = f(x);
    splines[0].x = x;
    for (int i = 1; i < N; i++, x += h) {
        splines[i].a = f(x);
        splines[i].c = c[i];
        splines[i].d = (c[i] - c[i - 1])/h;
        splines[i].b = h*c[i]/2 - h*h*splines[i].d/6 + (splines[i].a - splines[i - 1].a)/h;
        splines[i].x = x;
    }
}

double spline(Spline spline, double x)
{
    return spline.a + spline.b*(x - spline.x) + spline.c/2*(x - spline.x)*(x - spline.x) +\
            spline.d/6*(x - spline.x)*(x - spline.x)*(x - spline.x);
}

int main()
{
    double matrix[N][N + 1], vector[N];
    Spline splines[N];

    get_matrix(matrix, A, B);
    Tridiagonal(matrix, vector);
    get_splines(A, B, vector, splines);

    std::ofstream file(FILE_PATH, std::ofstream::out);

    int i = 0;
    const double h = (double)(B - A)/N;
    double y;
    for (double x = A; x <= B; x += h, i++) {
        y = spline(splines[i], x);
        file << x << ";" << y << std::endl;
        printf("%.4f, %f\n", x, y);
    }

    file.close();

    printf("\nFile \"%s\" has been successfully generated.\n", FILE_PATH);

    return 0;
}
