#include "stdmultilevelsolver.h"

namespace Gascoigne
{
  
  class MyMLS : public StdMultiLevelSolver
  {
    
  public:
    double NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo);
    void newton(Matrix& A, VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info);
    
    ///--------------- NUR Fuer FV transprt 
    
    void AssembleMatrix(Matrix& A, VectorInterface& u)
    {
      GetSolver()->MatrixZero(A);
      GetSolver()->AssembleMatrix(A,u,1.);
    }
    void ComputeIlu(Matrix& A, VectorInterface &u)
    {
      GetSolver()->ComputeIlu(A, u);
    }
    ///--------------- NUR Fuer FV transprt 
 
    
  };
}
