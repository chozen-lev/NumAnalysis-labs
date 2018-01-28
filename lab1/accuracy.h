#ifndef ACCURACY_H
#define ACCURACY_H

namespace accuracy {
    struct Result
    {
        double eps;
        int n;
        double value;
        double R;
    };

    Result exp(const double x, const double eps = 1e-3, int n = 1);
}

#endif // ACCURACY_H
