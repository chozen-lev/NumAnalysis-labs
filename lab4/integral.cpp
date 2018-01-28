#include "integral.h"
#include <cmath>

Result integral(double a, double b, func_t f, double h, int n)
{
    double In = 0.0;

    for (int i = 1; i < n; i++) {
        In += f(a + i*h);
    }

    In += (f(a) + f(a + n*h))/2;

    return Result { h, n, In*h };
}

Result RefinedCalculation(double a, double b, func_t f, double eps)
{
    double In = 0.0, I2n = 0.0, h = sqrt(eps);
    int n = (b - a)/h + 1;
    h = (b - a)/n;

    In = integral(a, b, f, h, n).value;
    n *= 2;
    h /= 2;
    I2n = integral(a, b, f, h, n).value;

    while (fabs(In - I2n) > 3*eps) {
        In = I2n;
        n *= 2;
        h /= 2;
        I2n = integral(a, b, f, h, n).value;
    }

    return Result { h, n, I2n };
}

double integralNL(double a, double b, func_t F)
{
    return F(b) - F(a);
}
