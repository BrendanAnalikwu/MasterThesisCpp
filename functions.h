#ifndef BENCHMARK_FUNCTIONS_H
#define BENCHMARK_FUNCTIONS_H

#include <cmath>
#include <map>
#include "paramfile.h"
#include "gascoigne.h"
#include "dataformathandler.h"
#include "filescanner.h"

namespace Gascoigne
{
class FourierSum
{
private:
    DoubleVector coef_x_c, coef_y_c, coef_x_s, coef_y_s;
    DoubleVector indices_x_c, indices_y_c, indices_x_s, indices_y_s;
public:

    std::map<double, double> map_x_c, map_y_c, map_x_s, map_y_s;
    FourierSum(const std::string& name, const ParamFile& pf)
    {
        DataFormatHandler DFH;
        DFH.insert(name + "_i_x_c", &indices_x_c);
        DFH.insert(name + "_x_c", &coef_x_c);
        DFH.insert(name + "_i_y_c", &indices_y_c);
        DFH.insert(name + "_y_c", &coef_y_c);
        DFH.insert(name + "_i_x_s", &indices_x_s);
        DFH.insert(name + "_x_s", &coef_x_s);
        DFH.insert(name + "_i_y_s", &indices_y_s);
        DFH.insert(name + "_y_s", &coef_y_s);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(pf, "Equation");
        assert(indices_x_c.size() == coef_x_c.size());
        assert(indices_y_c.size() == coef_y_c.size());
        assert(indices_x_s.size() == coef_x_s.size());
        assert(indices_y_s.size() == coef_y_s.size());

        fill_mapping(map_x_c, indices_x_c, coef_x_c);
        fill_mapping(map_y_c, indices_y_c, coef_y_c);
        fill_mapping(map_x_s, indices_x_s, coef_x_s);
        fill_mapping(map_y_s, indices_y_s, coef_y_s);
    }

    double inline operator()(double x, double y) const
    {
        return fourier_sum(map_x_c, map_x_s, x) * fourier_sum(map_y_c, map_y_s, y);
    }

    static void fill_mapping(std::map<double, double>& mapping, const DoubleVector& indices, const DoubleVector& coef)
    {
        for (int i = 0; i < indices.size(); ++i) { mapping[indices[i]] = coef[i]; }
    }

    static double fourier_sum(const std::map<double, double>& coef_cos, const std::map<double, double>& coef_sin, double x)
    {
        double res = 0.;
        for (auto it: coef_cos) { res += it.second * cos(it.first * x); }
        for (auto it: coef_sin) { res += it.second * sin(it.first * x); }
        return res;
    }

};
}

#endif //BENCHMARK_FUNCTIONS_H
