/*
 *
 * Lab2, "РОЗВ’ЯЗАННЯ РІВНЯНЬ З ОДНИМ НЕВІДОМИМ"
 * Lev Bazylskyi, KV-51
 * 22. 2*sqrt(2+3*sin(x))-x/2=0, (І,Д)
*/

#include <iostream>
#include <math.h>

double f(double x)
{
    return 2*sqrt(2+3*sin(x))-x/2;
}

double df(double x)
{
    return 3*cos(x)/sqrt(2+3*sin(x))-0.5;
}

double phi(double x, double lambda)
{
    return df(x) > 0.0 ? x - lambda*f(x) : x + lambda*f(x);
}

double iterativeMethod(double a, double b, double eps, double &delta)
{
    double m1 = fabs(df(a)), M1 = fabs(df(b));
    if (m1 > M1) {
        std::swap(m1, M1);
    }
    double x0, xk = (a+b)/2, lambda = 1.0/M1, q = 1.0 - m1/M1;

    do
    {
        x0 = xk;
        xk = phi(x0, lambda);
    } while (fabs(xk-x0) > (1-q)/q*eps);
    delta = fabs(xk-x0);

    return xk;
}

double tangentsMethod(double a, double b, double eps, double &delta)
{
    double m1 = fabs(df(a)), M1 = fabs(df(b));
    if (m1 > M1) {
        std::swap(m1, M1);
    }
    double xk = a;
    do
    {
        xk = xk - f(xk)/df(xk);
    } while (fabs(f(xk))/m1 > eps);
    delta = fabs(f(xk))/m1;

    return xk;
}

void printTableIterative(double a, double b)
{
    double x, delta = 0.0;
    printf("[%.1f;%.1f], iterative\neps\tX\t\t\tdelta\n", a, b);
    for (double eps = 1e-2; eps > 1e-14; eps *= 1e-3)
    {
        x = iterativeMethod(a, b, eps, delta);
        printf("%.0e\t%.14f\t%.6e\n", eps, x, delta);
    }
}

void printTableTangents(double a, double b)
{
    double x, delta = 0.0;
    printf("[%.1f;%.1f], tangents\neps\tX\t\t\tdelta\n", a, b);
    for (double eps = 1e-2; eps > 1e-14; eps *= 1e-3)
    {
        x = tangentsMethod(a, b, eps, delta);
        printf("%.0e\t%.14f\t%.6e\n", eps, x, delta);
    }
}

int main()
{
    printTableIterative(3.2, 3.7); printf("\n");
    printTableIterative(6.2, 6.7); printf("\n");
    printTableIterative(8.2, 8.7); printf("\n");

    printTableTangents(3.2, 3.7); printf("\n");
    printTableTangents(6.2, 6.7); printf("\n");
    printTableTangents(8.2, 8.7); printf("\n");

    return 0;
}
