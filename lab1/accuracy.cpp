#include "accuracy.h"
#define _USE_MATH_DEFINES 
#include <math.h>

double u(double x, int k)
{
    if (k >= 1) {
        return x / k * u(x, k - 1);
    }

    return 1.0;
}

double sum(double x, int k)
{
    if (k >= 1) {
        return sum(x, k - 1) + u(x, k);
    }

    return 1.0;
}

accuracy::Result accuracy::exp(const double x, const double eps, int n)
{
    double intpart, fractpart, R, exponential, value = 1.0;

    fractpart = modf(x, &intpart);
    exponential = x > 0.0 ? M_E : 1.0 / M_E;

    for (int i = 0; i < fabs(intpart); ++i) {
        value *= exponential;
    }

    if (eps != 0.0) {
        while (fabs(u(fractpart, n)) > eps) {
            n++;
        }
    }

    value *= sum(fractpart, n);
    R = u(fractpart, n);

    return Result { eps, n, value, R };
}
