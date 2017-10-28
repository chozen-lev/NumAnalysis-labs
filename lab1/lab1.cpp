/*
 *
 * Lab1, "ОБЧИСЛЕННЯ ЗНАЧЕНЬ ФУНКЦІЇ"
 * Lev Bazylskyi, KV-51
 * 22. e^x [-9.3;15]
*/

#include <iostream>
#include <math.h>

#define a (-9.3)
#define b 15.0
#define h (b-a)/10

double exponential(int x)
{
    double result = 1.0, exp = x > 0 ? M_E : 1/M_E;
    for (int i = 0; i < abs(x); i++) {
        result *= exp;
    }
    return result;
}

double get_u(int k, double x)
{
    return k >= 1 ? x/k*get_u(k-1, x) : 1.0;
}

double get_sum(int k, double x)
{
    return k >= 1 ? get_sum(k-1, x)+get_u(k, x) : 1.0;
}

int get_n(double eps, double x)
{
    int k = 0;
    do
    {
        k++;
    } while(abs(get_u(k, x)) >= eps);
    return k;
}

// int fact(int n)
// {
//     return n > 1 ? n*fact(n - 1) : 1;
// }

//double get_r(int n, double x)
//{
//    return pow(abs(x), n+1)/fact(n)/n;
//}

int main()
{
    double x = (b+a)/2, intpart, fractpart = modf(x, &intpart);
    int n, N;

    printf("eps\tn\tdelta\t\tR\n");
    for (double eps = 1e-2; eps >= 1e-14; eps *= 1e-3)
    {
        n = get_n(eps, fractpart);
        printf("%.0e\t%d\t%.6e\t%.6e\n", eps, n, abs(exponential(int(intpart))*get_sum(n, fractpart) - exp(x)), get_u(n,/*x*/ fractpart));

        if (eps == 1e-8) {
            N = n;
        }
    }

    printf("\n x\t\tdelta\t\tR\n");
    for (int i = 0; i <= 10; i++)
    {
        x = a+h*i;
        fractpart = modf(x, &intpart);
        printf("%c%.6f\t%.6e\t%.6e\n", x >= 0 ? ' ' : '\0', x, abs(exponential(int(intpart))*get_sum(N, fractpart) - exp(x)), get_u(N,/*x*/ fractpart));
    }

    return 0;
}
