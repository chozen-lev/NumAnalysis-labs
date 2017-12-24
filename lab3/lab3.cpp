/*
 *
 * Lab3, "РОЗВ’ЯЗАННЯ СИСТЕМ ЛІНІЙНИХ АЛГЕБРАЇЧНИХ РІВНЯНЬ"
 * Lev Bazylskyi, KV-51
 * 22. (22 % 6 = 4). Виключення Гаусса-Жордана,  метод ітерації Зейделя
*/

#include <iostream>
#include <math.h>
#include <cstring>

#define N 4
#define EPS 0.001

#define MIN -5
#define MAX 5

double const vector[N] = { 128.0, 325.0, 240.0, 120.0 };
double const matrix[N][N] = { { 10.0, 11.0, 14.0, 3.0 },
                              { 12.0, 45.0, 17.0, 15.0 },
                              { 9.0, 16.0, 31.0, 5.0 },
                              { 0.0, 4.0, 17.0, 15.0 } };

void print_matrix(const double matrix[][N], const double vector[N])
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%0.1f\t", matrix[i][j]);
        }
        printf("%0.1f\n", vector[i]);
    }
}

void print_answer(const double vector[N])
{
    for (int i = 0; i < N; i++) {
        printf("x%d = %.1f\n", i + 1, vector[i]);
    }
}

void GJ(double matrix[][N], double vector[N])
{
    double buff;

    for (int k = 0; k < N; k++) {
        if (matrix[k][k] == 0) {
            continue;
        }

        buff = matrix[k][k];
        for (int l = 0; l < N; l++) {
            matrix[k][l] /= buff;
        }
        vector[k] /= buff;

        for (int i = 0; i < N; i++) {
            if (i == k) {
                continue;
            }

            buff = matrix[i][k];
            for (int l = 0; l < N; l++) {
                matrix[i][l] -= matrix[k][l] * buff;
            }
            vector[i] -= vector[k] * buff;
        }
    }
}

int Seidel(double matrix[][N], double vector[N], double eps)
{
    int iter = 0;
    double q = 0.0, buff, norm;
    double alfa[N][N], beta[N], Xn[N];

    for (int i = 0; i < N; i++) {
        beta[i] = vector[i]/matrix[i][i];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                alfa[i][j] = -matrix[i][j]/matrix[i][i];
            } else {
                alfa[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        buff = 0;
        for (int j = 0; j < N; j++) {
            buff += fabs(alfa[i][j]);
        }
        if (buff > q) {
            q = buff;
        }
    }
//    printf("\nALFA-BETA check:\n"); print_matrix(alfa, beta);

    do {
        iter++;
//        printf("%d.\n", iter); print_answer(vector); printf("\n");

        for (int i = 0; i < N; i++) {
            Xn[i] = beta[i];
            for (int j = 0; j < i; j++) {
                Xn[i] += alfa[i][j]*Xn[j];
            }
            for (int j = i; j < N; j++) {
                Xn[i] += alfa[i][j]*vector[j];
            }
        }

        norm = 0.0;

        for (int i = 0; i < N; i++) {
            norm += (Xn[i] - vector[i])*(Xn[i] - vector[i]);
        }
        norm = sqrt(norm);

        if (norm <= eps*(1 - q)/q) {
            break;
        }

        memcpy(vector, Xn, N*sizeof(double));
    } while (true);

    return iter;
}

int main()
{
    double buff[N][N], buff2[N];

    print_matrix(matrix, vector);

    memcpy(buff, matrix, N * N * sizeof(double));
    memcpy(buff2, vector, N * sizeof(double));
    GJ(buff, buff2);
    printf("\nGauss-Jordan\n");
    print_answer(buff2);

//    double a1, a2, a3, a4, b;
//    for (int i = MIN; i < MAX; i++) {
//        for (int ii = MIN; ii < MAX; ii++) {
//            for (int j = MIN; j < MAX; j++) {
//                for (int jj = MIN; jj < MAX; jj++) {
//                    a1 = matrix[0][0] * i + matrix[1][0] * ii + matrix[2][0] * j + matrix[3][0] * jj;
//                    a2 = matrix[0][1] * i + matrix[1][1] * ii + matrix[2][1] * j + matrix[3][1] * jj;
//                    a3 = matrix[0][2] * i + matrix[1][2] * ii + matrix[2][2] * j + matrix[3][2] * jj;
//                    a4 = matrix[0][3] * i + matrix[1][3] * ii + matrix[2][3] * j + matrix[3][3] * jj;
//                    if (fabs(a4) > fabs(a2) + fabs(a3) + fabs(a1)) {
//                        b = vector[0] * i + vector[1] * ii + vector[2] * j + vector[3] * jj;
//                        printf("%d %d %d %d\n%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n\n", i, ii, j, jj, a1, a2, a3, a4, b);
//                    }
//                }
//            }
//        }
//    }

    memcpy(buff, matrix, N * N * sizeof(double));
    memcpy(buff2, vector, N * sizeof(double));

    // (1) -> 2*(1) - (3)
    buff[0][0] = 2*matrix[0][0] - matrix[2][0];
    buff[0][1] = 2*matrix[0][1] - matrix[2][1];
    buff[0][2] = 2*matrix[0][2] - matrix[2][2];
    buff[0][3] = 2*matrix[0][3] - matrix[2][3];
    buff2[0] = 2*vector[0] - vector[2];
    // (4) -> (3) - 2*(4)
    buff[3][0] = matrix[2][0] - 2*matrix[3][0];
    buff[3][1] = matrix[2][1] - 2*matrix[3][1];
    buff[3][2] = matrix[2][2] - 2*matrix[3][2];
    buff[3][3] = matrix[2][3] - 2*matrix[3][3];
    buff2[3] = vector[2] - 2*vector[3];
//    printf("\nTransformed matrix:\n"); print_matrix(buff, buff2);

    int iter = Seidel(buff, buff2, EPS);
    printf("\nSeidel, eps = %.0e\n", EPS);
    print_answer(buff2);
    printf("Iterations = %d\n", iter);

    return 0;
}
