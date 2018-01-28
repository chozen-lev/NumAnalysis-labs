#include "equation.h"
#include <math.h>

Result iterative(double a, double b, func_t f, func_t df, double eps)
{
    double m1 = fabs(df(a)), M1 = fabs(df(b));
    if (m1 > M1) {
        double buff = m1;
        m1 = M1;
        M1 = buff;
    }
    double x0, xk = (a + b)/2.0, lambda = 1.0/M1, q = 1.0 - m1/M1, delta;
    int iter = 0;

    do
    {
        x0 = xk;
        xk = df(x0) > 0.0 ? x0 - lambda*f(x0) : x0 + lambda*f(x0);
        iter++;
    } while (fabs(xk - x0) > (1 - q)/q*eps);
    delta = fabs(xk - x0)*q/(1 - q);

    return Result { xk, delta, iter };
}

Result tangents(double a, double b, func_t f, func_t df, func_t d2f, double eps)
{
    double m1 = fmin(fabs(df(a)), fabs(df(b))), xk, delta;
    xk = f(a)*d2f(a) > 0 ? a : b;
    int iter = 0;

    do
    {
        iter++;
        xk = xk - f(xk)/df(xk);
    } while (fabs(f(xk))/m1 > eps);
    delta = fabs(f(xk))/m1;

    return Result { xk, delta, iter };
}
