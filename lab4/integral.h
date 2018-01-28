#ifndef INTEGRAL_H
#define INTEGRAL_H

typedef double (*func_t)(double);

typedef struct
{
    double step;
    int n;
    double value;
} Result;

Result integral(double a, double b, func_t f, double h, int n);
Result RefinedCalculation(double a, double b, func_t f, double eps = 1e-3);
double integralNL(double a, double b, func_t F);

#endif // INTEGRAL_H
