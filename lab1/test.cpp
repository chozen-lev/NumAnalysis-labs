#include "accuracy.h"
#include <iostream>
#include <math.h>

const double A = -9.3;
const double B = 15.0;

const double eps_start = 1e-2;
const double eps_end = 1e-14;
const double eps_step = 1e-3;

int print_table1(const double eps_res = 1e-8)
{
    accuracy::Result result;
    double x = (B + A) / 2.0;
    int n = -1;

    printf("eps\tn\tdelta\t\tR\n");
    for (double eps = eps_start; eps >= eps_end; eps *= eps_step) {
        result = accuracy::exp(x, eps);
        if (eps == eps_res) {
            n = result.n;
        }

        printf("%.0e\t%d\t%.6e\t%.6e\n", result.eps, result.n, fabs(result.value - exp(x)), result.R);
    }

    return n;
}

void print_table2(const int n)
{
    accuracy::Result result;
    double x = (B + A) / 2.0, h = (B - A)/10;

    printf("\nx\t\tdelta\t\tR\n");
    for (int i = 0; i <= 10; i++)
    {
        x = A + h*i;
        result = accuracy::exp(x, 0.0, n);
        printf("%.6f\t%.6e\t%.6e\n", x, fabs(result.value - exp(x)), result.R);
    }
}

int main()
{
    int n = print_table1();
    print_table2(n);

    return 0;
}
