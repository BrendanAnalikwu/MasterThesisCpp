/*----------------------------   equation.h     ---------------------------*/
/*      $Id: myequation.h,v 1.1 2009/10/26 08:42:26 richter Exp $                 */
#ifndef __equation_H
#define __equation_H
/*----------------------------   equation.h     ---------------------------*/


#define  __MyEquation_h

#include  "paramfile.h"
#include  "equation.h"
#include  "boundaryequation.h"
#include <fstream>

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{


class Field
{
public:
    int _NX, _NY;
    double _minx, _maxx, _miny, _maxy;
    vector<double> _DAT;

    void init(const string& fname, int NX, int NY,
              double minx, double miny, double maxx, double maxy)
    {
        _NX = NX;
        _NY = NY;
        _minx = minx;
        _miny = miny;
        _maxx = maxx;
        _maxy = maxy;
        ifstream in(fname.c_str(), std::ifstream::binary);
        assert(in.is_open());

        _DAT.resize(_NX * _NY);
        for (int x = 0; x < _NX; ++x)
            for (int y = 0; y < _NY; ++y)
                in >> _DAT[y * _NX + x];
        string tut;
        in >> tut;
        cout << tut << endl;
        assert(tut == "x");
        in.close();
    }

    double operator()(double x, double y) const
    {
        double X = (x - _minx) / (_maxx - _minx);
        double Y = (y - _miny) / (_maxy - _miny);
        assert(X >= 0);
        assert(X < 1);
        assert(Y >= 0);
        assert(Y < 1);
        int iX = static_cast<int> (X * _NX);
        int iY = static_cast<int> (Y * _NY);
        assert(iX >= 0);
        assert(iY >= 0);
        assert(iX < _NX);
        assert(iY < _NY);
        return _DAT[iY * _NX + iX];
    }
};


class SeaIceData
{
public:
    double alpha0, tau0, Cstern, Pstern, Tref, Lref, MZ, ellipse, C, Cdw, rho, rhow, theta_w, f;
    double deltamin, shock0;
    double vin0;

    // aus der rechten seite
    double rhoa, Cda, theta_a, windx, windy;

    void BasicInit(const ParamFile& pf);
};


// tut
class MyEquation : public virtual Equation
{

protected:
    SeaIceData data;

    mutable double shock;
    mutable double alpha, tau;
    mutable double visc;
    mutable double uwx;
    mutable double uwy;

    mutable FemFunction* oldu;
    mutable CellFunction* DGH;
    mutable double h_;

    mutable double v_in;
    mutable Vertex2d _n;


public:

    void SetFemData(FemData& q) const
    {
        assert(q.find("oldu") != q.end());
        oldu = &q["oldu"];

    }

    void SetCellData(CellData& q) const
    {
        assert(q.find("DGH") != q.end());
        DGH = &q["DGH"];
    }

    void point(double h, const FemFunction& U, const Vertex2d& v) const;

    MyEquation() { abort(); }

    MyEquation(const SeaIceData& dat) : data(dat) {}

    MyEquation* createNew() const { return new MyEquation(data); }

    std::string GetName() const { return "MyEquation"; }

    int GetNcomp() const { return 2; }

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
};

}


/*----------------------------   equation.h     ---------------------------*/
/* end of #ifndef __equation_H */
#endif
/*----------------------------   equation.h     ---------------------------*/
