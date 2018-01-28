#ifndef EQUATION_H
#define EQUATION_H

typedef double (*func_t)(double);

typedef struct
{
    double value;
    double delta;
    int n;
} Result;

typedef enum {
    ITERATIVE,
    TANGENTS
} Methods;

Result iterative(double a, double b, func_t f, func_t df, double eps);
Result tangents(double a, double b, func_t f, func_t df, func_t d2f, double eps);

#endif // EQUATION_H
