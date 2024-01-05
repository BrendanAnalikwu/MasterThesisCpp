#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "compose_name.h"
#include <algorithm>
#include "functions.h"

using namespace std;

#include "stopwatch.h"


extern ofstream ELLIPSE_OUT;

double TIME, DT, DTSUB, DELTAMIN, DTUL, starttime;
extern int zahl;
extern double STEUERUNG_MU;

namespace Gascoigne
{
extern Timer GlobalTimer;

void Loop::FVStep(const vector <vector<int>>& FV,
                  const vector <vector<Vertex2d>>& FV_midpoint,
                  GlobalVector& H,
                  const GlobalVector& V,
                  double dtFV, int M)
{
    GlobalVector HOLD = H;  // Copy H
    for (int iy = 0; iy < M; ++iy)
        for (int ix = 0; ix < M; ++ix)
        {
            int q = FV[ix][iy];
            const IntVector& ioc = GetMultiLevelSolver()->GetSolver()->GetMesh()->IndicesOfCell(q);

            for (int c = 0; c < 2; ++c)
            {
                H(q, c) = HOLD(q, c);

                //  2  3
                //  0  1

                if (ix > 0) // to the left
                {
                    double vx = .5 * V(ioc[0], 0) + .5 * V(ioc[2], 0);
                    int ql = FV[ix - 1][iy];
                    double flux = vx * ((vx < 0) ? HOLD(q, c) : HOLD(ql, c));
                    H(q, c) += flux * M / 0.5 * dtFV;
                }
                if (ix < M - 1) // to the right
                {
                    double vx = .5 * V(ioc[3], 0) + .5 * V(ioc[1], 0);
                    int qr = FV[ix + 1][iy];
                    double flux = - vx * ((vx > 0) ? HOLD(q, c) : HOLD(qr, c));
                    H(q, c) += flux * M / 0.5 * dtFV;
                }
                if (iy > 0) // down
                {
                    double vy = .5 * V(ioc[0], 1) + .5 * V(ioc[1], 1);
                    int qd = FV[ix][iy - 1];
                    double flux = vy * ((vy < 0) ? HOLD(q, c) : HOLD(qd, c));
                    H(q, c) += flux * M / 0.5 * dtFV;
                }
                if (iy < M - 1) // up
                {
                    double vy = .5 * V(ioc[2], 1) + .5 * V(ioc[3], 1);
                    int qu = FV[ix][iy + 1];
                    double flux = - vy * ((vy > 0) ? HOLD(q, c) : HOLD(qu, c));
                    H(q, c) += flux * M / 0.5 * dtFV;
                }
            }
        }

    // Restrict A to [0,1] and H >= 0
    for (int i = 0; i < H.n(); ++i)
    {
        H(i, 1) = min(1.0, H(i, 1));
        H(i, 1) = max(0.0, H(i, 1));
        H(i, 0) = max(0.0, H(i, 0));
    }
}


void Loop::run(const std::string& problemlabel)
{

    GlobalTimer.start("alles");
    GlobalTimer.start("-> Vorbereitung");


    zahl = 1;
    double tref;
    double dtmax, endtime;
    int prerefine;
    string _reloadu, _reloadh, _reloadoldu;
    string coef_name;
    bool _visudgh, _visuv0;

    // Load parameters from parameter file
    // Each block needs a different scope
    {
        DataFormatHandler DFH;
        DFH.insert("dt", &DT, 0.);
        DFH.insert("dtmax", &dtmax, 0.);
        DFH.insert("time", &starttime, 0.);
        DFH.insert("endtime", &endtime, 0.);
        DFH.insert("reloadu", &_reloadu);
        DFH.insert("reloadh", &_reloadh);
        DFH.insert("reloadoldu", &_reloadoldu);
        DFH.insert("visudgh", &_visudgh, false);
        DFH.insert("visuv0", &_visuv0, false);
        FileScanner FS(DFH);
        FS.NoComplain();

        FS.readfile(_paramfile, "Loop");
        assert(DT > 0.0);
        cout << starttime << "start" << endl;
    }
    {
        DataFormatHandler DFH;
        DFH.insert("prerefine", &prerefine, 0);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Mesh");
    }
    {
        DataFormatHandler DFH;
        DFH.insert("deltamin", &DELTAMIN, 0.);
        DFH.insert("Tref", &tref, 0.);
        DFH.insert("coef_name", &coef_name, "coef.param");
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Equation");
        assert(DELTAMIN > 0.0);
        assert(tref > 0.0);
    }
    {
        string discname;
        DataFormatHandler DFH;
        DFH.insert("discname", &discname);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Solver");
        assert(discname == "CGQ1" || discname == "CGQ2");
    }

    ParamFile coef_params(_s_resultsdir + string("/") + JOB_ARRAY_ID + string("/") + coef_name);
    FourierSum H_initial("H", coef_params, 0.), A_initial("A", coef_params, 0., 1.);

    // vectors for solution and right hand side
    Vector u("u"), f("f"), oldu("oldu"), other("other");
    Matrix A("A");

    // vector for ice height and concentration
    Vector DGH("DGH", "cell");

    PrintMeshInformation();
    // initialize problem, solver and vectors within Gascoigne
    GetMultiLevelSolver()->ReInit();
    GetMultiLevelSolver()->SetProblem("seaice");
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(oldu);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->ReInitVector(other);
    GetMultiLevelSolver()->ReInitVector(DGH);

    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    GetMultiLevelSolver()->GetSolver()->OutputSettings();

    // Finite Volume grid for advection
    const GascoigneMesh* M2 = dynamic_cast<const GascoigneMesh2d* >(GetMultiLevelSolver()->GetSolver()->GetMesh());
    assert(M2);

    int ncells = M2->ncells();  // Number of cells in grid
    int M = sqrt(ncells);  // Number of cells in each direction
    std::cout << "ncells: " << ncells << " " << M << std::endl;
    assert(ncells == M * M);

    // Initialise FV and FV_midpoint
    vector<vector<int> > FV(M, vector<int>(M)); // FV[i][j] is the cell number in the discretisation
    vector<vector<Vertex2d> > FV_midpoint(M, vector<Vertex2d>(M)); // FV_midpoint[i][j] is the midpoint of cel ij
    map<Vertex2d, int> MidPointToCellNumber;

    // Initialise H and A
    GlobalVector& glDGH = GetMultiLevelSolver()->GetSolver()->GetGV(DGH);

    // Fill mappings by looping over cells
    for (int q = 0; q < ncells; ++q)
    {
        // Get vertex indices of cell q
        const IntVector& ioc = M2->IndicesOfCell(q);  // {4*q, 1+4*q, 2+4*q, 3+4*q}
        // Average the vertex coordinates to compute cell midpoint
        Vertex2d v;
        v.equ(0.25, M2->vertex2d(ioc[0]), 0.25, M2->vertex2d(ioc[1]), 0.25,
              M2->vertex2d(ioc[2]), 0.25, M2->vertex2d(ioc[3]));

        // Fill mappings
        MidPointToCellNumber[v] = q;

        double x = v.x();
        double y = v.y();

        // fill H and A with initial values
        glDGH(q, 0) = H_initial(x, y);
        glDGH(q, 1) = A_initial(x, y);
    }

    // Fill FV and FV_midpoint mappings
    int ix = 0;
    int iy = 0;
    for (const auto& it: MidPointToCellNumber)
    {
        FV[ix][iy] = it.second;
        FV_midpoint[ix][iy] = it.first;

        // Increase ix, iy indices
        ++iy;
        if (iy == M)
        {
            iy = 0;
            ++ix;
        }
    }
    // End of finite volume setup


    // Momentum equation solving
    cout << "------------------------------" << endl;
    cout << "sea-ice " << endl;
    cout << "------------------------------" << endl << endl;

    // Initialise momentum with initial conditions
//    InitSolution(u);
    GlobalVector& gu = GetMultiLevelSolver()->GetSolver()->GetGV(u);
    FourierSum u0_1("u00", coef_params), u0_2("u01", coef_params);
    for (int i=0; i<M2->nnodes(); i++) {
        Vertex2d v = M2->vertex2d(i);
        gu(i, 0) = u0_1(v);
        gu(i, 1) = u0_2(v);
    }

    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
    GetMultiLevelSolver()->GetSolver()->SubtractMean(u);  // Don't know why this. TODO
    if (_visuv0) GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir + string("/") + JOB_ARRAY_ID + "/v0", u, 0);
    GetMultiLevelSolver()->Equ(oldu, 1.0, u);  // Save u to oldu


    double writenext = 0;
    string res;

    nvector<double> functionals;

    TIME = starttime;


    nvector<double> Jtotal, Jtotal1;

    GlobalTimer.stop("-> Vorbereitung");
    for (_iter = 1; _iter <= _niter; _iter++)
    {
        GlobalTimer.start("-> Iteration");

        // Increment time (dimensionless)
        TIME += DT;

        GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

        cout << "\n time step " << _iter << " "
             << "\t" << TIME * tref / 60 / 60 / 24 << " days" << endl;

        // FV-transport of sea ice thickness H and concentration A
//        GlobalTimer.start("--> Transport");
//        // Number of subcycles
//        int NSUB = 1;
//
//        double dtFV = DT / NSUB; // Finite volumes time step size
//        // Perform FV step
//        for (int ii = 1; ii <= NSUB; ++ii)
//            FVStep(FV, FV_midpoint, glDGH, GetMultiLevelSolver()->GetSolver()->GetGV(u), dtFV, M);
//
        if (_visudgh) GetMultiLevelSolver()->GetSolver()->CellVisu(_s_resultsdir + string("/") + JOB_ARRAY_ID + "/dgh", glDGH, _iter); // Save H & A to disc
//
//        GlobalTimer.stop("--> Transport");

        // End FV-Transport

        // Set sea ice problem (to solve sea ice momentum equation)
        cout << "Momentum" << endl;
        GlobalTimer.start("--> Momenten");
        GetMultiLevelSolver()->SetProblem("seaice");
        GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
        GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
        GetMultiLevelSolver()->GetSolver()->AddCellVector("DGH", DGH);
        GetMultiLevelSolver()->GetSolver()->Equ(oldu, 1.0, u); // Save previous value for u to oldu
        GetMultiLevelSolver()->AssembleMatrix(A, u);
        GetMultiLevelSolver()->ComputeIlu(A, u);

        res = Solve(A, u, f, "Results/u");
        assert(res == "converged");

        // Visualise sea ice velocity
        //  if(_iter==_niter){
        GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir + string("/") + JOB_ARRAY_ID + "/v", u, _iter); // Save new u to disc
        // }
        functionals = Functionals(u, f);

        GetMultiLevelSolver()->DeleteNodeVector("oldu");
        GetMultiLevelSolver()->GetSolver()->DeleteCellVector("DGH");
        GlobalTimer.stop("--> Momenten");

        writenext += DT;

        GlobalTimer.stop("-> Iteration");
        if (TIME >= endtime * 24 * 60 * 60 / tref) break;

    }

    GlobalTimer.stop("alles");

    GlobalTimer.print();
}

}






 
