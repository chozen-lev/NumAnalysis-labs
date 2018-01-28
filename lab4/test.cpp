#include "integral.h"
#include <iostream>
#include <cmath>
#include <vector>

const double A = 2.0;
const double B = 87.0;
const double M2 = 11.048309686793592; // x = 18.970216999999998

const double eps_start = 1e-2;
const double eps_end = 1e-8;
const double eps_step = 1e-2;

double f(double x) { return cos(x)*(1 + 2*exp(sin(x))); }
double F(double x) { return 2*exp(sin(x)) + sin(x); }

std::vector<double> print_table1()
{
    std::vector<double> deltas;
    Result result;
    double h, integ = integralNL(A, B, F);
    int n;

    printf("eps\th\tintegral\t\tdelta\n");
    for (double eps = eps_start; eps > eps_end; eps *= eps_step)
    {
        h = sqrt(12*eps/(B - A)/M2);
        n = (B - A)/h + 1;
        h = (B - A)/n;
        result = integral(A, B, f, h, n);
        deltas.push_back(fabs(integ - result.value));
        printf("%.0e\t%.1e\t%.14f\t%.4e\n", eps, result.step, result.value, deltas.back());
    }

    return deltas;
}

void print_table2(std::vector<double> deltas)
{
    Result result;
    double integ = integralNL(A, B, F);

    printf("\neps\t\th\tdelta\n");
    for (unsigned i = 0; i < deltas.size(); i++)
    {
        result = RefinedCalculation(A, B, f, deltas[i]);
        printf("%.4e\t%.1e\t%.4e\n", deltas[i], result.step, fabs(integ - result.value));
    }
}

int main()
{
    printf("integral = %.14f\n\n", integralNL(A, B, F));
    std::vector<double> deltas;
    deltas = print_table1();
    print_table2(deltas);

    return 0;
}
