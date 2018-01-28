#include <iostream>
#include "linear_system.h"

const double matrix[N][N + 1] = { { 10.0, 11.0, 14.0, 3.0, 128.0 },
                                  { 12.0, 45.0, 17.0, 15.0, 325.0 },
                                  { 9.0, 16.0, 31.0, 5.0, 240.0 },
                                  { 0.0, 4.0, 17.0, 15.0, 120.0 } };

void print_matrix(const double matrix[N][N + 1])
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N + 1; j++) {
            printf("%.1f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_roots(const double roots[N])
{
    for (int i = 0; i < N; ++i) {
        printf("x%d = %.1f\n", i + 1, roots[i]);
    }
}

void print_GaussJordan()
{
    double roots[N];
    GaussJordan(matrix, roots);
    printf("Gauss - Jordan method:\n");
    print_roots(roots);
    printf("\n");
}

void print_Seidel()
{
    double matrix2[N][N + 1];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N + 1; j++) {
            matrix2[i][j] = matrix[i][j];
        }
    }

    for (int i = 0; i < N + 1; ++i) {
        matrix2[0][i] = 2*matrix[0][i] - matrix[2][i];
        matrix2[3][i] = matrix[2][i] - 2*matrix[3][i];
    }

    double roots[N];
    Seidel(matrix2, roots);
    printf("Seidel method:\n");
    print_roots(roots);
    printf("\n");
}

int main()
{
    print_matrix(matrix);
    print_GaussJordan();
    print_Seidel();
    return 0;
}
