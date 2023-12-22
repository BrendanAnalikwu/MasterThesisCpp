#ifndef BENCHMARK_FUNCTIONS_H
#define BENCHMARK_FUNCTIONS_H

#include <cmath>
#include <map>

namespace Gascoigne
{
double fourier_sum(std::map<double, double> coef_cos, std::map<double, double> coef_sin, double x);
}


#endif //BENCHMARK_FUNCTIONS_H
