#/*----------------------------   problem.h     ---------------------------*/
/*      $Id: problem.h,v 1.2 2009/10/26 09:07:16 richter Exp $                 */
#ifndef __problem_H
#define __problem_H
/*----------------------------   problem.h     ---------------------------*/


#include  "problemdescriptorbase.h"
#include  "myequation.h"
#include  "dirichletdata.h"
#include  "domainrighthandside.h"
#include  "periodicdata.h"
#include<fstream>
#include  "zerodirichletdata.h"
#include "componentinformationbase.h"

extern double TIME, DT,DTSUB,DDD;

using namespace std;

namespace Gascoigne
{
  double WindX(double x, double y, double t)
  {
    //    return (x-250.)/20.0;
    
    double tP=t; while (tP>8) { tP-=8.0; }
    int p = ((static_cast<int>(t/8)))%2;
    double vmax = 15.0; // maximale windgeschwindigkeit in m/s
    double mx,my;
    double alpha = M_PI/2.0; // 90 grad ist ohne Konvergenz oder Divergenz
    // windstarke
    // double ws = tanh(tP*(8.0-tP)/2.0);  
    double ws = tanh(tP*(8.0-tP)/2.0);  

    if ((p%2) == 0)
      {
	mx = 50+800/16.0*tP;
	my = 50+800/16.0*tP;
	alpha -=  M_PI/2.0/5; // 18 Grad Konvergenz

      }
    else
      { ws=-ws;
	mx = 450-800/16.0*tP;
	my = 450-800/16.0*tP;
	alpha -=  M_PI/2.0/10; // 9 Grad Divergenz 10
	

      }

    double wx = cos(alpha)*(x-mx) + sin(alpha)*(y-my);
    //    double wy = -sin(alpha)*(x-mx) + cos(alpha)*(y-my);
    double r = sqrt((mx-x)*(mx-x)+(my-y)*(my-y));
    double s = 1.0/50.0*exp(-0.01*r);

      return -wx*s*ws*vmax;
  }
  double WindY(double x, double y, double t)
  {
    //return (y-250.)/20.0;

    double tP=t; while (tP>8) { tP-=8.0; }
    int p = ((static_cast<int>(t/8)))%2;
    double vmax = 15.0; // maximale windgeschwindigkeit
    double mx,my;
    double alpha = M_PI/2.0;  // 90 grad ist ohne Konvergenz oder Divergenz
    // windstarke

    //  double ws = 1.0;
    double ws=tanh(tP*(8.0-tP)/2.0);

    if ((p%2) == 0)
      {
	mx = 50+800/16.0*tP;
	my = 50+800/16.0*tP;
	alpha -=  M_PI/2.0/5; // 18 Grad Konvergenz
		

      }
    else
      {ws=-ws;
	mx = 450-800/16.0*tP;
	my = 450-800/16.0*tP;
	alpha -=  M_PI/2.0/10; // 9 Grad Divergenz
		

      }

    double wy = -sin(alpha)*(x-mx) + cos(alpha)*(y-my);
    double r = sqrt((mx-x)*(mx-x)+(my-y)*(my-y));
    double s = 1.0/50.0*exp(-0.01*r);


    
    
   return -wy * s * ws*vmax;
  }

  

  // class for specifying Dirichlet data
  class SeaIceDirichletData : public DirichletData
  {
  public:
    SeaIceDirichletData(const ParamFile& pf) : DirichletData(pf)
    {}
    std::string GetName() const { return "SeaIce Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
    }
  };

  

  
  /*---------------------------------------------------*/





  // class for specifying problem's right hand side
  class MyRightHandSide : public DomainRightHandSide
  {
    SeaIceData data;
    
  public:

   
    void SetFemData(FemData& q) const 
    {
    }

    int GetNcomp()        const { return 4; }
    std::string GetName() const { return "Rhs"; }
    

  MyRightHandSide(const SeaIceData& dat) : DomainRightHandSide(), data(dat)
      {
	assert(data.rho>0);
	assert(data.rhoa>0);
	assert(data.Lref>0);
	assert(data.Tref>0);

      }

     MyRightHandSide* createNew() const { return new MyRightHandSide(data); }

    double operator()(int comp, const Vertex2d& v) const
    {
     
      
      double SZ = data.rhoa/data.rho * data.Cda * data.Lref;
      
      double time = TIME * data.Tref;
      double X = v.x()*data.Lref;
      double Y = v.y()*data.Lref;
    
      
      double U_x = WindX(X/1000, Y/1000, TIME*data.Tref/(60.0*60.0*24.0)); // umrechnen in km 
      double U_y = WindY(X/1000, Y/1000, TIME*data.Tref/(60.0*60.0*24.0));// umrechnen in km 
 

      double ux = U_x/data.Lref * data.Tref;
      double uy = U_y/data.Lref * data.Tref;
      
 
      double vw    = sqrt(ux*ux + uy*uy);
    

     
      if (comp==0)
	return SZ * vw * (ux * cos(data.theta_a) - uy * sin(data.theta_a));
      else if (comp==1)
	return SZ * vw * (uy * cos(data.theta_a) + ux * sin(data.theta_a));
      else
	return 0;
   
    }
  };
  

  /*--------------------------------------------------*/
  
  
  class MyCI : public ComponentInformationBase 
  {
    
  public:
      
    virtual void      GetScalarName   (IndexType i, std::string& s_name)
    {
      if (i==0) s_name = "vx";
      if (i==1) s_name = "vy";
    }
    virtual const IndexType GetNVectors     () const { return 1; }
    virtual void      GetVectorName   (IndexType i, std::string& s_name) const { s_name = "V";}
    virtual void      GetVectorIndices(IndexType i, std::array<int,3>& fa) const
    {
      fa[0]=0; fa[1]=1; fa[2]=-1;
    }
  };

  
  class MyCIWater : public ComponentInformationBase 
  {
    
  public:

    virtual void      GetScalarName   (IndexType i, std::string& s_name)
    {
      if (i==0) s_name = "p";
      if (i==1) s_name = "vx";
      if (i==2) s_name = "vy";
    }
    virtual const IndexType GetNVectors     () const { return 1; }
    virtual void      GetVectorName   (IndexType i, std::string& s_name) const { s_name = "V";}
    virtual void      GetVectorIndices(IndexType i, std::array<int,3>& fa) const
    {
      fa[0]=1; fa[1]=2; fa[2]=-1;
    }
  };

  
  class MyCITr : public ComponentInformationBase 
  {
    
  public:

    virtual void      GetScalarName   (IndexType i, std::string& s_name)
    {
      if (i==0) s_name = "h";
      if (i==1) s_name = "A";
    }
    virtual const IndexType GetNVectors() const { return 0; }
  };



  // main class for defining the problem to solve
  class SeaIceProblem : public ProblemDescriptorBase
  {
    SeaIceData data;
  public:
    Application* NewRightHandSide() const
    {
      return new MyRightHandSide(data);
    }

    Equation *NewEquation() const
    {
      return new MyEquation(data);
    }
    
    
  
    std::string GetName() const {return "SeaIce Problem";}
    void BasicInit(const ParamFile& pf)
    {
      data.BasicInit(pf);
      
      // equation to solve
      GetEquationPointer()      = new MyEquation(data);

      // definition of dirichlet boundary data
      //	GetDirichletDataPointer() = new ZeroDirichletData;
      GetDirichletDataPointer() = new SeaIceDirichletData(pf);
      GetComponentInformationPointer() = new MyCI;
      GetRightHandSidePointer()=new MyRightHandSide(data);
      // 
      ProblemDescriptorBase::BasicInit(pf);
    }
  };
  



  
  // class for specifying Dirichlet data
  class OtherDirichletData : public DirichletData
  {
  public:
    OtherDirichletData(const ParamFile& pf) : DirichletData(pf)
    {}

    std::string GetName() const { return "Transport Dirichlet Data"; }
    
    void operator()(DoubleVector& b, const Vertex2d& v, int col) const
    {
      b.zero();
    }

  };

  // main class for defining the problem to solve
  class TGProblem : public ProblemDescriptorBase
  {
  public:
    SeaIceData data;
    
    std::string GetName() const {return "TG Problem";}

    
    void BasicInit(const ParamFile& pf)
    {
      data.BasicInit(pf);
     
      GetEquationPointer()      = new TransportEquation(pf);// Not in use !!!!!
	
      ProblemDescriptorBase::BasicInit(pf);

    }
  };

  // to visualize the sheat stress

  class OtherProblem : public ProblemDescriptorBase
  {
    SeaIceData data;
  public:

    Equation *NewEquation() const
    {
      return new OtherEquation(data);
    }

  
    std::string GetName() const {return "Other Problem";}
    void BasicInit(const ParamFile& pf)
    {
      data.BasicInit(pf);
      
      // equation to solve
      GetEquationPointer()      = new OtherEquation(data);
      GetDirichletDataPointer() = new OtherDirichletData(pf);

      ProblemDescriptorBase::BasicInit(pf);
      //GetPeriodicDataPointer() = new MyPeriodic;
    }
  };


  
  class MyBM : public BoundaryManager
  {
    IntSet _nocolors;
  public:
    const IntSet& GetDirichletDataColors        () const
    { return _nocolors; }
  };

}

/*----------------------------   problem.h     ---------------------------*/
/* end of #ifndef __problem_H */
#endif
/*----------------------------   problem.h     ---------------------------*/
  
