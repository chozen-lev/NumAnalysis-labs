#include "equation.h"
#include <iostream>
#include <math.h>
#include <vector>

const double eps_start = 1e-2;
const double eps_end = 1e-14;
const double eps_step = 1e-3;

double f(double x) { return 2*sqrt(2 + 3*sin(x)) - x/2; }
double df(double x) { return 3*cos(x)/sqrt(2 + 3*sin(x)) - 0.5; }
double d2f(double x) { return (-3*sin(x)*sqrt(2+3*sin(x))-4.5*cos(x)*cos(x)/sqrt(2+3*sin(x)))/(2+3*sin(x)); }

std::vector<Result> print_table(double a, double b, Methods method)
{
    std::vector<Result> res;
    printf("[%.2f;%.2f], %s\neps\tX\t\t\tdelta\n", a, b, method == ITERATIVE ? "iterative" : "tangents");
    for (double eps = eps_start; eps > eps_end; eps *= eps_step)
    {
        switch (method) {
        case ITERATIVE:
            res.push_back(iterative(a, b, f, df, eps));
            break;
        case TANGENTS:
            res.push_back(tangents(a, b, f, df, d2f, eps));
            break;
        }
        printf("%.0e\t%.14f\t%.6e\n", eps, res.back().value, res.back().delta);
    }
    printf("\n");

    return res;
}

void print_compare(double a, double b, std::vector<Result> iterative, std::vector<Result> tangents)
{
    double eps = eps_start;
    printf("[%.1f;%.1f], compare\neps\tIterative\tTangents\n", a, b);
    for (unsigned i = 0; i < iterative.size(); ++i, eps *= eps_step) {
        printf("%.0e\t%d\t\t%d\n", eps, iterative[i].n, tangents[i].n);
    }
    printf("\n");
}

int main()
{
    std::vector<Result> iterative, tangents;
    iterative = print_table(3.2, 3.7, ITERATIVE);
    tangents = print_table(3.2, 3.7, TANGENTS);
    print_table(6.2, 6.7, ITERATIVE);
    print_table(6.2, 6.7, TANGENTS);
    print_table(8.2, 8.7, ITERATIVE);
    print_table(8.2, 8.7, TANGENTS);
    print_compare(3.2, 3.7, iterative, tangents);

    return 0;
}
