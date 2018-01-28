#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#define N 4

void GaussJordan(const double matrix[][N + 1], double roots[N]);
void Seidel(const double matrix[][N + 1], double roots[N], double eps = 1e-3);

#endif // LINEAR_SYSTEM_H
