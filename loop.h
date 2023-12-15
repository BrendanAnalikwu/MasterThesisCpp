/*----------------------------   loop.h     ---------------------------*/
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/


#include "stdloop.h"
#include  "periodicmapping.h"
#include "meshagent.h"
#include "gascoignemesh2d.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "simplematrix.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include "solver.h"


namespace Gascoigne
{


class Loop : public StdLoop
{

public:
    void BasicInit(const ParamFile& paramfile, const ProblemContainer* PC, const FunctionalContainer* FC)
    {


        GetMultiLevelSolverPointer() = new MyMLS();

        StdLoop::BasicInit(paramfile, PC, FC);
    }

    void FVStep(const std::vector<std::vector<int> >& FV,
                const std::vector<std::vector<Vertex2d> >& FV_midpoint,
                GlobalVector& H,
                const GlobalVector& V,
                double dtFV, int M);

    std::string PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s);

    void run(const std::string& problemlabel);


};

}



/*----------------------------   loop.h     ---------------------------*/
#endif
/*----------------------------   loop.h     ---------------------------*/
