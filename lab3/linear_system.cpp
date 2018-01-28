#include "linear_system.h"
#include <cstring>
#include <cmath>

void GaussJordan(const double _matrix[][N + 1], double roots[N])
{
    double matrix[N][N + 1], buff;

    for (int i = 0; i < N; ++i) {
        memcpy(matrix[i], _matrix[i], (N + 1)*sizeof(double));
    }

    for (int k = 0; k < N; k++) {
        if (matrix[k][k] == 0) {
            continue;
        }

        buff = matrix[k][k];
        for (int l = 0; l < N + 1; l++) {
            matrix[k][l] /= buff;
        }

        for (int i = 0; i < N; ++i) {
            if (i == k) {
                continue;
            }

            buff = matrix[i][k];
            for (int l = 0; l < N + 1; l++) {
                matrix[i][l] -= matrix[k][l] * buff;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        roots[i] = matrix[i][N];
    }
}

void Seidel(const double matrix[][N + 1], double roots[N], double eps)
{
    double alfa[N][N];
    double beta[N], Xprev[N];
    double q = 0.0, norm = 0.0;

    for (int i = 0; i < N; ++i) {
        beta[i] = matrix[i][N]/matrix[i][i];
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                alfa[i][j] = -matrix[i][j]/matrix[i][i];
            } else {
                alfa[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += fabs(alfa[i][j]);
        }
        if (sum > q) {
            q = sum;
        }
    }

    for (int i = 0; i < N; ++i) {
        Xprev[i] = beta[i];
        for (int j = 0; j < i; j++) {
            Xprev[i] += alfa[i][j]*Xprev[j];
        }
        for (int j = i; j < N; j++) {
            Xprev[i] += alfa[i][j]*matrix[j][N];
        }
    }

    while (true) {
        for (int i = 0; i < N; ++i) {
            roots[i] = beta[i];
            for (int j = 0; j < i; j++) {
                roots[i] += alfa[i][j]*roots[j];
            }
            for (int j = i; j < N; j++) {
                roots[i] += alfa[i][j]*Xprev[j];
            }
        }

        norm = 0.0;

        for (int i = 0; i < N; ++i) {
            norm += (roots[i] - Xprev[i])*(roots[i] - Xprev[i]);
        }
        norm = sqrt(norm);

        if (norm <= eps*(1 - q)/q) {
            break;
        }

        for (int i = 0; i < N; ++i) {
            Xprev[i] = roots[i];
        }
    }
}
