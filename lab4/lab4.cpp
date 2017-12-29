/*
 *
 * Lab4, "ЧИСЕЛЬНЕ ІНТЕГРУВАННЯ"
 * Lev Bazylskyi, KV-51
 * 22. integrate cos(x)*(1+2*e^sin(x)) from 2 to 87, F = 2*e^sin(x)+sin(x)
 * (22 % 2 = 0). Узагальнена формула трапецій
*/

#include <iostream>
#include <math.h>
#include <vector>

#define A 2
#define B 87

#define EPS_START 1e-2
#define EPS_END 1e-8
#define EPS_STEP 1e-2

//#define d2f(x) (2*cos(x)*(exp(sin(x))*cos(x)*cos(x)-exp(sin(x))-(2*exp(sin(x)+1)*cos(x)-4*exp(sin(x))*sin(x)*cos(x))))
//#define M2 abs(d2f(18.970216999999998))
#define M2 11.048309686793592 // x = 18.970216999999998

double f(const double x)
{
    return cos(x)*(1 + 2*exp(sin(x)));
}

double F(const double x)
{
    return 2*exp(sin(x)) + sin(x);
}

double Trapezoidal(const double a, const double b, double h)
{
    double In = 0.0;
    int n = (b - a)/h;
    h = (b - a)/n;

    for (int i = 1; i < n; i++) {
        In += f(a + i*h);
    }

    In += (f(a) + f(b))/2;

    return In*h;
}

double print_table1(const double integral, const double eps)
{
    double value, delta, h = sqrt(12*eps/(B - A)/M2);

    value = Trapezoidal(A, B, h);
    delta = abs(integral - value);

    printf("%.0e\t%.1e\t%.12f\t%.4e\n", eps, h, value, delta);

    return delta;
}

void print_table2(const double integral, const double eps)
{
    long long int n = (B - A)/sqrt(eps);
    double In = 0.0, I2n = 0.0;

    In = Trapezoidal(A, B, (double)(B - A)/n); // h = (B - A)/n
    n <<= 1;
    I2n = Trapezoidal(A, B, (double)(B - A)/n);

    while (fabs(In - I2n) > 3*eps) {
        In = I2n;
        n <<= 1;
        I2n = Trapezoidal(A, B, (double)(B - A)/n);
    }

    In = I2n;

    printf("%.4e\t%.1e\t%.4e\n", eps, (double)(B - A)/n, abs(integral - In));
}

int main()
{
    std::vector<double> deltas;
    double IntegralNL = F(B) - F(A);
    printf("Newton–Leibniz= %.12f\n\n", IntegralNL);

    printf("eps\th\tintegral\tdelta\n");
    for (double eps = EPS_START; eps >= EPS_END; eps *= EPS_STEP) {
        deltas.push_back(print_table1(IntegralNL, eps));
    }

    printf("\neps\t\th\tdelta\n");
    for (unsigned i = 0; i < deltas.size(); i++) {
        print_table2(IntegralNL, deltas[i]);
    }

    return 0;
}
