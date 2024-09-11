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
    double scaling_min{}, scaling_max{}, scaling{}, intercept{}, lb{}, ub{};
    bool lb_set{}, ub_set{};
public:

    std::map<double, double> map_x_c, map_y_c, map_x_s, map_y_s;

    FourierSum() = default;

    FourierSum(const std::string& name, const ParamFile& pf, const double lb, const double ub) : lb(lb), ub(ub)
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
        DFH.insert(name + "_min", &scaling_min);
        DFH.insert(name + "_max", &scaling_max);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(pf, "Coefficients");
        assert(indices_x_c.size() == coef_x_c.size());
        assert(indices_y_c.size() == coef_y_c.size());
        assert(indices_x_s.size() == coef_x_s.size());
        assert(indices_y_s.size() == coef_y_s.size());

        double sigma_x = 0;
        for (const auto& c: coef_x_c) sigma_x += std::abs(c);
        for (const auto& c: coef_x_s) sigma_x += std::abs(c);
        double sigma_y = 0;
        for (const auto& c: coef_y_c) sigma_y += std::abs(c);
        for (const auto& c: coef_y_s) sigma_y += std::abs(c);
        double sigma = sigma_x * sigma_y;

        double range = scaling_max - scaling_min;
        intercept = range / 2. + scaling_min;
        scaling = range / (2. * sigma);
//        std::cout << scaling << ", " << intercept  << ", " << coef_x_c[0] << std::endl;

        fill_mapping(map_x_c, indices_x_c, coef_x_c);
        fill_mapping(map_y_c, indices_y_c, coef_y_c);
        fill_mapping(map_x_s, indices_x_s, coef_x_s);
        fill_mapping(map_y_s, indices_y_s, coef_y_s);

        lb_set = true;
        ub_set = true;
    }

    FourierSum(const std::string& name, const ParamFile& pf, const double lb): FourierSum(name, pf, lb, 0.) {
        ub_set = false;
    }

    FourierSum(const std::string& name, const ParamFile& pf): FourierSum(name, pf, 0., 0.) {
        lb_set = false;
        ub_set = false;
    }

    double inline operator()(double x, double y) const
    {
        double res = fourier_sum(map_x_c, map_x_s, x) * fourier_sum(map_y_c, map_y_s, y);
        res *= scaling;
        res += intercept;
        if(lb_set) res = std::max(lb, res);
        if(ub_set) res = std::min(ub, res);
        return res;
    }

    double inline operator()(const Vertex2d& v) const { return (*this)(v.x(), v.y()); }

    static void fill_mapping(std::map<double, double>& mapping, const DoubleVector& indices, const DoubleVector& coef)
    {
        for (int i = 0; i < indices.size(); ++i) { mapping[indices[i]] = coef[i]; }
    }

    static double inline
    fourier_sum(const std::map<double, double>& coef_cos, const std::map<double, double>& coef_sin, double x)
    {
        double res = 0.;
        for (auto it: coef_cos) { res += it.second * cos(it.first * x); }
        for (auto it: coef_sin) { res += it.second * sin(it.first * x); }
        return res;
    }

};

class Cyclone
{
private:
    double c;
public:
    double mx, my, ws, vmax, alpha, r0, vx_m, vy_m;

    Cyclone() = default;

    Cyclone(const ParamFile& pf)
    {
        bool cyclone;
        DataFormatHandler DFH;
        DFH.insert("W_mx", &mx, .256);
        DFH.insert("W_my", &my, .256);
        DFH.insert("W_vx_m", &vx_m, .0005787037);  // 4e-4 - 8e-4
        DFH.insert("W_vy_m", &vy_m, .0005787037);
        DFH.insert("W_vmax", &vmax, .015);
        DFH.insert("W_alpha", &alpha, M_PI * 2. / 5.);
        DFH.insert("W_r0", &r0, .1);
        DFH.insert("W_cyclone", &cyclone);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(pf, "Cyclone");

        if (cyclone) ws = 1.;
        else ws = -1.;
        c = vmax * ws / r0 / exp(-1.);
    }

    double inline x(const Vertex2d& v, const double t) const { return x(v.x(), v.y(), t); }

    double inline y(const Vertex2d& v, const double t) const { return y(v.x(), v.y(), t); }

    double inline x(const double x, const double y, const double t) const
    {
        double mx_t = mx + vx_m * t;
        double my_t = my + vy_m * t;

        return c * (cos(alpha) * (x - mx_t) + sin(alpha) * (y - my_t)) *
           exp(-sqrt((x - mx_t) * (x - mx_t) + (y - my_t) * (y - my_t)) / r0);
    }

    double inline y(const double x, const double y, const double t) const
    {
        double mx_t = mx + vx_m * t;
        double my_t = my + vy_m * t;
        return c * (-sin(alpha) * (x - mx_t) + cos(alpha) * (y - my_t)) *
               exp(-sqrt((x - mx_t) * (x - mx_t) + (y - my_t) * (y - my_t)) / r0);
    }
};
}

#endif //BENCHMARK_FUNCTIONS_H
