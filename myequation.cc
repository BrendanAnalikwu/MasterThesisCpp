#include  "myequation.h"
#include  "filescanner.h"


extern double TIME;
extern double DT, DDD;
extern double DELTAMIN;
extern double STEUERUNG_MU, zahl;


/*-----------------------------------------*/

namespace Gascoigne
{
void SeaIceData::BasicInit(const ParamFile& paramfile)
{
    string coef_name, _s_resultsdir;

    {
        DataFormatHandler DFH;
        DFH.insert("v_in", &vin0, 0.);
        DFH.insert("rho", &rho, 0.);
        DFH.insert("rhow", &rhow, 0.);
        DFH.insert("Tref", &Tref, 0.0);
        DFH.insert("Lref", &Lref, 0.0);
        DFH.insert("Pstern", &Pstern, 2.75e4);
        DFH.insert("ellipse", &ellipse, 2.0);
        DFH.insert("C", &C, 20.0);
        DFH.insert("Cdw", &Cdw, 5.2e-3);
        DFH.insert("f", &f, 0.0);
        DFH.insert("theta_w", &theta_w, 0.0);

        DFH.insert("rhoa", &rhoa, 0.);          // aus Rhs
        DFH.insert("Cda", &Cda, 0.0);
        DFH.insert("theta_a", &theta_a, 0.0);
        DFH.insert("windx", &windx, 0.0);
        DFH.insert("windy", &windy, 0.0);

        DFH.insert("coef_name", &coef_name, "coef.param");

        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(paramfile, "Equation");
    }
    {
        DataFormatHandler DFH;
        DFH.insert("resultsdir", &_s_resultsdir, "Results");

        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(paramfile, "Loop");
    }
    assert(rho > 0);
    assert(Tref > 0);
    assert(Lref > 0);
    assert(Pstern > 0);
    assert(ellipse > 0);
    assert(C > 0);
    //assert(f>0);

    MZ = 0.5 * Tref * Tref * Pstern / rho / Lref / Lref;
    std::cout << "Mehlmann-Zahl " << MZ << std::endl;

    ParamFile coef_params(_s_resultsdir + string("/") + JOB_ARRAY_ID + string("/") + coef_name);
    Windx = FourierSum("Wx", coef_params);
    Windy = FourierSum("Wy", coef_params);
    Oceanx = FourierSum("Ox", coef_params);
    Oceany = FourierSum("Oy", coef_params);
}

// sea ice mometum equation
void MyEquation::point(double h, const FemFunction& U, const Vertex2d& v) const
{
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    h_ = h;

    uwx = data.Oceanx(v);
    uwy = data.Oceany(v);
}

/*-----------------------------------------*/

// u[0], u[1] v
//tutu
// u[2] h    u[3] a

void MyEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
    // ocean drag tau_w
    // wind drag is caculated in problem.h 
    double WZ = data.rhow / data.rho * data.Cdw * data.Lref;

    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;

    double vd = sqrt(vd_x * vd_x + vd_y * vd_y + 1.e-8);

    b[0] += WZ * vd * (vd_x * cos(data.theta_w) - vd_y * sin(data.theta_w)) * N.m();
    b[1] += WZ * vd * (vd_y * cos(data.theta_w) + vd_x * sin(data.theta_w)) * N.m();


    double dmin = DELTAMIN * data.Tref;

    // DELTA of the viscosity 
    double DELTAsquare =
            (1.0 + pow(data.ellipse, -2.0)) * (U[0].x() * U[0].x() + U[1].y() * U[1].y())
            + pow(data.ellipse, -2.0) * pow(U[0].y() + U[1].x(), 2.0)
            + 2.0 * (1.0 - pow(data.ellipse, -2.0)) * U[0].x() * U[1].y();

    DELTAsquare += dmin * dmin;
    //normal 
    double DELTA = sqrt(DELTAsquare);


    //---------------FV----------------------------------------- 
    // term Correolis 1.46 e -4 /s
    b[0] += -data.Tref * (*DGH)[0] * data.f * (vd_y) * N.m();
    b[1] += data.Tref * (*DGH)[0] * data.f * (vd_x) * N.m();


    // time derivative
    b[0] += (*DGH)[0] * (U[0].m() - (*oldu)[0].m()) / DT * N.m();
    b[1] += (*DGH)[0] * (U[1].m() - (*oldu)[1].m()) / DT * N.m();


    // part of the ice strength P
    double ef = exp(-data.C * (1.0 - (*DGH)[1]));

    // ice strength P
    b[0] += -data.MZ * (*DGH)[0] * ef * N.x();
    b[1] += -data.MZ * (*DGH)[0] * ef * N.y();

    /////  (sigma, nabla phi)	
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
            b[i] += data.MZ * (*DGH)[0] * pow(data.ellipse, -2.0) * ef / DELTA * (U[i][j + 1] + U[j][i + 1]) * N[j + 1];
            b[i] += data.MZ * (*DGH)[0] * (1.0 - pow(data.ellipse, -2.0)) * ef / DELTA * U[j][j + 1] * N[i + 1];
        }

    //---------------FV-----------------------------------------


}
/*-----------------------------------------*/
// sea ice mometum equation Matrix
void MyEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{

    //---------------------------------------
    //  tau_w ocean drag 
    double WZ = data.rhow / data.rho * data.Cdw * data.Lref;

    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;

    double vd = sqrt(vd_x * vd_x + vd_y * vd_y + 1.e-8);
    double vd_0 = 1.0 / (2.0 * vd) * 2.0 * vd_x * M.m();
    double vd_1 = 1.0 / (2.0 * vd) * 2.0 * vd_y * M.m();

    //b[0] += WZ * vd   *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0, 0) += WZ * vd_0 * (vd_x * cos(data.theta_w) - vd_y * sin(data.theta_w)) * N.m();
    A(0, 1) += WZ * vd_1 * (vd_x * cos(data.theta_w) - vd_y * sin(data.theta_w)) * N.m();
    A(0, 0) += WZ * vd * (M.m() * cos(data.theta_w)) * N.m();
    A(0, 1) += WZ * vd * (-M.m() * sin(data.theta_w)) * N.m();

    //b[1] += WZ * vd   * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1, 0) += WZ * vd_0 * (vd_y * cos(data.theta_w) + vd_x * sin(data.theta_w)) * N.m();
    A(1, 1) += WZ * vd_1 * (vd_y * cos(data.theta_w) + vd_x * sin(data.theta_w)) * N.m();
    A(1, 1) += WZ * vd * (M.m() * cos(data.theta_w)) * N.m();
    A(1, 0) += WZ * vd * (M.m() * sin(data.theta_w)) * N.m();

    //---------------------------------------

    //---------------------------------------
    //derivative of Delta
    double dmin = DELTAMIN * data.Tref;
    double DELTAsquare =
            (1.0 + pow(data.ellipse, -2.0)) * (U[0].x() * U[0].x() + U[1].y() * U[1].y())
            + pow(data.ellipse, -2.0) * pow(U[0].y() + U[1].x(), 2.0)
            + 2.0 * (1.0 - pow(data.ellipse, -2.0)) * U[0].x() * U[1].y();
    DELTAsquare += dmin * dmin;
    double DELTA = sqrt(DELTAsquare);

    // derivativ w.r.t Delta
    double DELTAsquare_0 =
            (1.0 + pow(data.ellipse, -2.0)) * (2.0 * U[0].x() * M.x())
            + pow(data.ellipse, -2.0) * 2.0 * (U[0].y() + U[1].x()) * M.y()
            + 2.0 * (1.0 - pow(data.ellipse, -2.0)) * M.x() * U[1].y();
    double DELTAsquare_1 =
            (1.0 + pow(data.ellipse, -2.0)) * (2.0 * U[1].y() * M.y())
            + pow(data.ellipse, -2.0) * 2.0 * (U[0].y() + U[1].x()) * M.x()
            + 2.0 * (1.0 - pow(data.ellipse, -2.0)) * U[0].x() * M.y();

    // For DDD=1.0 one gets the standard Newotn method 
    double DELTA_0 = -0.5 * pow(DELTAsquare, -1.5) * DELTAsquare_0 * DDD;
    double DELTA_1 = -0.5 * pow(DELTAsquare, -1.5) * DELTAsquare_1 * DDD;
    //---------------------------------------


    //--------------------------------------
    //Derivative Coriolis term 
    //b[0] += -data.Tref*(*H)[0].m()*f*(vd_y)*N.m();
    A(0, 1) += -data.Tref * (*DGH)[0] * data.f * M.m() * N.m();
    //b[1] +=  data.Tref*(*H)[0].m()*f*(vd_x)* N.m();
    A(1, 0) += data.Tref * (*DGH)[0] * data.f * M.m() * N.m();
    //--------------------------------------

    //--------------------------------------
    // derivative time dependen term 
    A(0, 0) += (*DGH)[0] * M.m() / DT * N.m();
    A(1, 1) += (*DGH)[0] * M.m() / DT * N.m();
    //--------------------------------------

    //--------------------------------------
    // derivative of the stress tensor

    double ef = exp(-data.C * (1.0 - (*DGH)[1]));

    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
            //	  b[i] += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
            A(i, 0) += data.MZ * (*DGH)[0] * pow(data.ellipse, -2.0) * ef * DELTA_0 * (U[i][j + 1] + U[j][i + 1]) *
                       N[j + 1];
            A(i, 1) += data.MZ * (*DGH)[0] * pow(data.ellipse, -2.0) * ef * DELTA_1 * (U[i][j + 1] + U[j][i + 1]) *
                       N[j + 1];

            A(i, i) += data.MZ * (*DGH)[0] * pow(data.ellipse, -2.0) * ef / DELTA * M[j + 1] * N[j + 1];
            A(i, j) += data.MZ * (*DGH)[0] * pow(data.ellipse, -2.0) * ef / DELTA * M[i + 1] * N[j + 1];


            A(i, 0) += data.MZ * (*DGH)[0] * (1.0 - pow(data.ellipse, -2.0)) * ef * DELTA_0 * U[j][j + 1] * N[i + 1];
            A(i, 1) += data.MZ * (*DGH)[0] * (1.0 - pow(data.ellipse, -2.0)) * ef * DELTA_1 * U[j][j + 1] * N[i + 1];

            // // b[i] += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N[i+1];
            A(i, j) += data.MZ * (*DGH)[0] * (1.0 - pow(data.ellipse, -2.0)) * ef / DELTA * M[j + 1] * N[i + 1];
        }

    // sigma_1
    // A(i,i) += data.MZ * (*DGH)[0]* pow(data.ellipse,-2.0) * ef / DELTA * M[j+1] * N[j+1];
    // A(i,j) += data.MZ * (*DGH)[0]* pow(data.ellipse,-2.0) * ef / DELTA * M[i+1] * N[j+1]
    // A(i,j) += data.MZ * (*DGH)[0]* (1.0-pow(data.ellipse,-2.0)) * ef / DELTA * M[j+1]*N[i+1];

    // sigma_2 derivative w.r.t. Delta
    //A(i,0) += data.MZ * (*DGH)[0]* pow(data.ellipse,-2.0) * ef * DELTA_0 * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
    // A(i,1) += data.MZ * (*DGH)[0]* pow(data.ellipse,-2.0) * ef * DELTA_1 * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
    //A(i,0) += data.MZ * (*DGH)[0] * (1.0-pow(data.ellipse,-2.0)) * ef * DELTA_0 * U[j][j+1]*N[i+1];
    //A(i,1) += data.MZ * (*DGH)[0]* (1.0-pow(data.ellipse,-2.0)) * ef * DELTA_1 * U[j][j+1]*N[i+1];
    //--------------------------------------



}

}
