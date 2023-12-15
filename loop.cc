#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "compose_name.h"
#include <algorithm>

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


    GlobalVector HOLD = H;
    for (int iy = 0; iy < M; ++iy)
        for (int ix = 0; ix < M; ++ix)
        {
            int q = FV[ix][iy];
            const IntVector& ioc = GetMultiLevelSolver()->GetSolver()->GetMesh()->IndicesOfCell(q);

            for (int c = 0; c < 2; ++c)
            {
                H(q, c) = HOLD(q, c);


                //
                //  2  3
                //
                //  0  1


                if (ix > 0) // nach links
                {
                    double vx = .5 * V(ioc[0], 0) + .5 * V(ioc[2], 0);
                    double vy = .5 * V(ioc[0], 1) + .5 * V(ioc[2], 1);
                    int ql = FV[ix - 1][iy];
                    double flux = HOLD(q, c) * (-1.0) * vx;
                    if (-1.0 * vx < 0)
                        flux = HOLD(ql, c) * (-1.0) * vx;

                    H(q, c) -= flux * M / 0.5 * dtFV;
                }
                if (ix < M - 1) // nach rechts
                {
                    double vx = .5 * V(ioc[3], 0) + .5 * V(ioc[1], 0);
                    double vy = .5 * V(ioc[3], 1) + .5 * V(ioc[1], 1);
                    int qr = FV[ix + 1][iy];
                    double flux = HOLD(q, c) * (1.0) * vx;
                    if (+1.0 * vx < 0)
                        flux = HOLD(qr, c) * (1.0) * vx;
                    H(q, c) -= flux * M / 0.5 * dtFV;
                }
                if (iy > 0) // nach unten
                {
                    double vx = .5 * V(ioc[0], 0) + .5 * V(ioc[1], 0);
                    double vy = .5 * V(ioc[0], 1) + .5 * V(ioc[1], 1);
                    int qu = FV[ix][iy - 1];
                    double flux = HOLD(q, c) * (-1.0) * vy;
                    if (-1.0 * vy < 0)
                        flux = HOLD(qu, c) * (-1.0) * vy;

                    H(q, c) -= flux * M / 0.5 * dtFV;
                }
                if (iy < M - 1) // nach oben
                {
                    double vx = .5 * V(ioc[2], 0) + .5 * V(ioc[3], 0);
                    double vy = .5 * V(ioc[2], 1) + .5 * V(ioc[3], 1);
                    int qo = FV[ix][iy + 1];
                    double flux = HOLD(q, c) * (1.0) * vy;
                    if (+1.0 * vy < 0)
                        flux = HOLD(qo, c) * (1.0) * vy;
                    H(q, c) -= flux * M / 0.5 * dtFV;
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
    if (1)
    {
        DataFormatHandler DFH;
        DFH.insert("dt", &DT, 0.);
        DFH.insert("dtmax", &dtmax, 0.);
        DFH.insert("time", &starttime, 0.);
        DFH.insert("endtime", &endtime, 0.);
        DFH.insert("reloadu", &_reloadu);
        DFH.insert("reloadh", &_reloadh);
        DFH.insert("reloadoldu", &_reloadoldu);
        FileScanner FS(DFH);
        FS.NoComplain();

        FS.readfile(_paramfile, "Loop");
        assert(DT > 0.0);
        cout << starttime << "start" << endl;
    }
    if (1)
    {
        DataFormatHandler DFH;
        DFH.insert("prerefine", &prerefine, 0);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Mesh");
    }
    if (1)
    {
        DataFormatHandler DFH;
        DFH.insert("deltamin", &DELTAMIN, 0.);
        DFH.insert("Tref", &tref, 0.);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Equation");
        assert(DELTAMIN > 0.0);
        assert(tref > 0.0);
    }
    string discname;
    if (1)
    {
        DataFormatHandler DFH;
        DFH.insert("discname", &discname);
        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(_paramfile, "Solver");
        assert(discname == "CGQ1" || discname == "CGQ2");
    }




    // vectors for solution and right hand side
    Vector u("u"), f("f"), oldu("oldu"), other("other");
    Matrix A("A");

    Vector DGH("DGH", "cell");

    for (int ADAITER = 0; ADAITER < 1; ++ADAITER)
    {


        PrintMeshInformation();
        // initialize problem, solver and vectors
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

        ///////////// Finite Volumen grid for advection
        const GascoigneMesh* M2 = dynamic_cast<const GascoigneMesh2d* >(GetMultiLevelSolver()->GetSolver()->GetMesh());
        assert(M2);

        int ncells = M2->ncells();  // Number of cells in grid
        int M = sqrt(ncells);  // Number of cells in each direction
        std::cout << "ncells: " << ncells << " " << M << std::endl;
        assert(ncells = M * M);

        // Initialise FV and FV_midpoint
        vector<vector<int> > FV(M, vector<int>(
                M)); // FV[i][j] is the cell number in the CR-discretisation. FV is a vector of vectors
        vector<vector<Vertex2d> > FV_midpoint(M, vector<Vertex2d>(M)); // FV[i][j] is the midpoint of cel ij

        // Initialise and fill midpoint to cell number mapping
        map<Vertex2d, int> MidPointToCellNumber;
        for (int q = 0; q < ncells; ++q)
        {
            const IntVector& ioc = M2->IndicesOfCell(q);  // {4*q, 1+4*q, 2+4*q, 3+4*q}
            Vertex2d v;
            v.equ(0.25, M2->vertex2d(ioc[0]), 0.25, M2->vertex2d(ioc[1]), 0.25, M2->vertex2d(ioc[2]), 0.25,
                  M2->vertex2d(ioc[3]));
            MidPointToCellNumber[v] = q;
        }
        int ix = 0;
        int iy = 0;
        for (auto it: MidPointToCellNumber)
        {
            FV[ix][iy] = it.second;
            ++iy;
            if (iy == M)
            {
                iy = 0;
                ++ix; ///////////// Finite Volumen Gitter
                const GascoigneMesh* M2 = dynamic_cast<const GascoigneMesh2d* >(GetMultiLevelSolver()->GetSolver()->GetMesh());
                assert(M2);

                int ncells = M2->ncells();
                int M = sqrt(ncells);  // anzahl zellen in jede Richtung
                assert(ncells = M * M);
                vector<vector<int> > FV(M, vector<int>(M)); // FV[i][j] ist nummer der Zelle in CR-Diskretisierung
                vector<vector<Vertex2d> > FV_midpoint(M, vector<Vertex2d>(M)); // FV[i][j] ist mittelpunkt der Zelle
                map<Vertex2d, int> MidPointToCellNumber;
                for (int q = 0; q < ncells; ++q)
                {
                    const IntVector& ioc = M2->IndicesOfCell(q);
                    Vertex2d v;

                    // Compute the midpoint of the cell by taking a linear combination of the grid points
                    v.equ(0.25, M2->vertex2d(ioc[0]), 0.25, M2->vertex2d(ioc[1]), 0.25, M2->vertex2d(ioc[2]), 0.25,
                          M2->vertex2d(ioc[3]));
                    MidPointToCellNumber[v] = q;
                }
                // Loop over each cell and fill FV (which is a map from i,j grid node numbers to the cell number
                int ix = 0;
                int iy = 0;
                for (auto it: MidPointToCellNumber)
                {
                    FV[ix][iy] = it.second;
                    ++iy;
                    if (iy == M)
                    {
                        iy = 0;
                        ++ix;
                    }
                }

                // Fill mapping from i,j grid node numbers to midpoint
                for (int iy = 0; iy < M; ++iy)
                    for (int ix = 0; ix < M; ++ix)
                    {
                        const IntVector& ioc = M2->IndicesOfCell(FV[ix][iy]);
                        Vertex2d v;
                        v.equ(0.25, M2->vertex2d(ioc[0]), 0.25, M2->vertex2d(ioc[1]), 0.25, M2->vertex2d(ioc[2]), 0.25,
                              M2->vertex2d(ioc[3]));
                        FV_midpoint[ix][iy] = v;
                    }

            }
        }



        // Inizialize H and A
        GlobalVector& glDGH = GetMultiLevelSolver()->GetSolver()->GetGV(DGH);
        for (int iy = 0; iy < M; ++iy)
            for (int ix = 0; ix < M; ++ix)
            {
                double x = FV_midpoint[ix][iy].x();
                double y = FV_midpoint[ix][iy].y();

                glDGH(FV[ix][iy], 0) = 0.005 * sin(60 * x) + 0.005 * sin(30 * y) + 0.3;
                glDGH(FV[ix][iy], 1) = 1.0;
            }

        //////////////////////////////////////////////////



        cout << "------------------------------" << endl;
        cout << "sea-ice " << endl;
        cout << "------------------------------" << endl << endl;

        // See-Eis
        GetMultiLevelSolver()->SetProblem("seaice");
        if (_initial == "reload")
        {
            abort();
            GetMultiLevelSolver()->GetSolver()->Read(u, _reloadu);
        } else
        {
            InitSolution(u);
            GetMultiLevelSolver()->Equ(oldu, 1.0, u);
        }


        double stepback = 0.0;
        double writenext = 0;
        int writeiter = 0;
        string res;

        nvector<double> functionals;


        int timeinc = 0;

        clock_t start, end;
        double cpu_time_used, a;

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




            //-----------------------------FV--transport of sea ice thickness H and concentration A
            GlobalTimer.start("--> Transport");
            // Number of subcycles
            int NSUB = 1;

            double dtFV = DT / NSUB; // Finite volumes time step size
            // Perform FV step
            for (int ii = 1; ii <= NSUB; ++ii)
                FVStep(FV, FV_midpoint, glDGH, GetMultiLevelSolver()->GetSolver()->GetGV(u), dtFV, M);

            GetMultiLevelSolver()->GetSolver()->CellVisu("Results/dgh", glDGH,
                                                         _iter + ADAITER * 1000); // Save H & A to disc

            GlobalTimer.stop("--> Transport");

            //-----------------------------FV--Transport----------------------------

            // set sea ice problem (to solve sea ice momentum equation)
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

            // Visulaize sea ice velocity
            //  if(_iter==_niter){
            GetMultiLevelSolver()->GetSolver()->Visu("Results/v", u, _iter + ADAITER * 1000); // Save new u to disc
            // }
            functionals = Functionals(u, f);

            GetMultiLevelSolver()->DeleteNodeVector("oldu");
            GetMultiLevelSolver()->GetSolver()->DeleteCellVector("DGH");
            GlobalTimer.stop("--> Momenten");

            GlobalTimer.start("--> Other");
            // Other-Problem ---> to visulize shear stress
            GetMultiLevelSolver()->SetProblem("other");
            GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
            GetMultiLevelSolver()->AddNodeVector("U", u);
            GetMultiLevelSolver()->GetSolver()->AddCellVector("DGH", DGH);
            Solve(A, other, f, "Results/o");

            //visulaize shear stress
            // if(_iter==_niter){
            GetMultiLevelSolver()->GetSolver()->Visu("Results/o", other, _iter + ADAITER * 1000);
            // }
            GetMultiLevelSolver()->DeleteNodeVector("U");
            GetMultiLevelSolver()->GetSolver()->DeleteCellVector("DGH");

            stringstream str;
            str << "fv_bench_kinetic.txt";
            ofstream OUTF(str.str().c_str(), ios::app);
            OUTF.precision(10);
            OUTF << TIME << " " << functionals << endl;
            GlobalTimer.stop("--> Other");


            writenext += DT;


            GlobalTimer.stop("-> Iteration");
            if (TIME >= endtime * 24 * 60 * 60 / tref) break;

        }

    }

    GlobalTimer.stop("alles");

    GlobalTimer.print();
}


}






 
