#include "functions.h"

namespace Gascoigne
{
double fourier_sum(std::map<double, double> coef_cos, std::map<double, double> coef_sin, double x)
{
    double res = 0.;

    for (auto it: coef_cos) { res += it.second * cos(it.first * x); }
    for (auto it: coef_sin) { res += it.second * sin(it.first * x); }

    return res;
}
}