#include "stdmultilevelsolver.h"

namespace Gascoigne
{

class MyMLS : public StdMultiLevelSolver
{

public:
    double NewtonUpdate(double& rr, Vector& x, Vector& dx, Vector& r, const Vector& f, NLInfo& nlinfo);

    void newton(Matrix& A, Vector& u, const Vector& f, Vector& r, Vector& w, NLInfo& info);

    ///--------------- NUR Fuer FV transprt 

    void AssembleMatrix(Matrix& A, Vector& u)
    {
        GetSolver()->MatrixZero(A);
        GetSolver()->AssembleMatrix(A, u, 1.);
    }

    void ComputeIlu(Matrix& A, Vector& u)
    {
        GetSolver()->ComputeIlu(A, u);
    }
    ///--------------- NUR Fuer FV transprt 


};
}
