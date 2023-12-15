
#include "solver.h"

int zahl, it;
double DDD = 1.0;
int writeiter;


using namespace std;

namespace Gascoigne
{
double MyMLS::NewtonUpdate(double& rr, Vector& x, Vector& dx, Vector& r, const Vector& f, NLInfo& nlinfo)
{
    const CGInfo& linfo = nlinfo.GetLinearInfo();
    bool lex = linfo.control().status() == "exploded";

    double nn = NewtonNorm(dx);

    double nnl2 = GetSolver()->GetGV(dx).norm();

    double nr = GetSolver(ComputeLevel)->Norm(r);

    ostringstream ss;
    ss << "mat_nn_it_DD_neu1" << it << ".dat";
    string datname = ss.str();
    /*
    ofstream OUT(datname.c_str());
    if(it==12){
     OUT << "Update Norm Lunendlich "<< nn
      << "\t update norm l2 " << nnl2
      << "\t norm residuum " << nr << endl;
     OUT.close();

    }
    */
    if (nn > 1.e30) lex = 1;
    if (!(nn >= 0.)) lex = 1;
    if (nr > 1.e30) lex = 1;
    if (!(nr >= 0.)) lex = 1;

    if (lex)
    {
        nlinfo.control().status() = "diverged";
        cerr << "linear : " << linfo.control().status() << endl;
        cerr << "nonlinear : " << nn << endl;
        return NewtonNorm(dx);
    }

    double omega = 0.25;
    double relax = 1.0;


    GetSolver(ComputeLevel)->SetPeriodicVectorZero(dx);
    /*
    // if(zahl==1000 && it==3)
    //{
    // nur test
    GlobalVector  X = GetSolver(ComputeLevel)->GetGV(x);
    NewtonResidual(r,x,f);
    double rold = NewtonNorm(r);

    for (double om = 0.0010406;om>0.0001;om*=0.99)
      {
        GetSolver(ComputeLevel)->GetGV(x)=X;
        GetSolver(ComputeLevel)->Add(x,om,dx);
        NewtonResidual(r,x,f);
        rr = NewtonNorm(r);
        cerr << om << "\t" << rr/rold << endl;
      }

    GetSolver(ComputeLevel)->GetGV(x)=X;
    //}
    */
    // So gehoerts:
    GetSolver(ComputeLevel)->Add(x, relax, dx);
    NewtonResidual(r, x, f);
    rr = NewtonNorm(r);

    string message = "";
    int diter = 0;


    for (diter = 0; diter < nlinfo.user().maxrelax(); diter++)
        // for (diter=0;diter<15;diter++)
    {
        message = nlinfo.check_damping(diter, rr);


        if (message == "ok") break;
        if (message == "continue")
        {
            GetSolver(ComputeLevel)->Add(x, relax * (omega - 1.), dx);

            NewtonResidual(r, x, f);
            rr = NewtonNorm(r);
            relax *= omega;
            if (it > 10)
            {
                //	  cout<<relax<<"relax_netwon"<<endl
                ;
            }
            continue;
        }
        if (message == "exploded")
        {
            GetSolver(ComputeLevel)->Add(x, -relax, dx);
            relax = 0.;
            cout << "Damping exploded !!!!!" << endl;
            nlinfo.control().status() = "diverged";
            break;
        }
    }

    nlinfo.check_damping(diter, rr);



    //  NewtonUpdateShowCompResiduals(nlinfo.control().iteration(), x, r, f,dx);

    return NewtonNorm(dx);
}


void MyMLS::newton(Matrix& A, Vector& u, const Vector& f, Vector& r, Vector& w, NLInfo& info)
{
    DDD = 1.0;

    double rho1 = 0.0;
    double rho2 = 0.0;
    double ET1 = 1.0;
    double ET2 = 1.0;
    info.reset();
    double rr = NewtonResidual(r, u, f);
    bool reached = info.check(0, rr, 0.);
    NewtonOutput(info);
    NewtonPreProcess(u, f, info);
    nvector<int> nlin, nresi;
    double res_old, res_new, res_mat;

    it = 0.0;

    std::vector<double> ressi(200);
    ressi[0] = 1.0;


    for (it = 1; !reached; it++)
    {
        NewtonMatrixControl(A, u, info);
        NewtonVectorZero(w);

        NewtonLinearSolve(A, w, r, info.GetLinearInfo());
        nlin.push_back(info.GetLinearInfo().control().iteration());

        double rw = NewtonUpdate(rr, u, w, r, f, info);

        ressi[it] = rr;
        reached = info.check(it, rr, rw);

        rho1 = info.statistics().lastrate();

        NewtonOutput(info);


        nresi.push_back(rr);

        ET1 = DDD;

        // DDD=1.0;
        //if(it==12)
        DDD = std::min(ET1 * (0.2 + 4 / (0.7 + exp(1.51 * rho1))), 1.);
        //cout<<it<<endl;
        //	cout<<"DDD"<<DDD<<"it"<<it<<endl;
        /*	//	cout<<"DDD"<<DDD<<"it"<<it<<endl;
            if (rho1>=1)
        {DDD=0.5;
          cout<<DDD<<endl;
        }
        if(fabs(rho1-1)<0.005)
          {	  DDD=1.0;
            cout<<DDD<<"DDD"<<endl;
          }
        //if (zahl==1 && it==2)
        //{DDD=0.5;
        //cout<<DDD<<"DDD"<<endl;
        // }
           //   DDD=max(DDD,0.1);

        */
        //	 if(DDD<0.3)

        if (DDD < 0.95)
            DDD = 1;
        cout << "Steuerung Newton rho/DDD/ET: " << rho1 << "\t" << DDD << "\t" << rr
             << endl; //	DDD=min(ET1*(0.2+4/(0.7+exp(1.51*rho1))),1);


    }


    // das gibt eine statistik ueber die Newton und ~GMRES konvergenz aus.
    //   cerr << "New: " << it <<"\t"<< nlin.sum() << "\t" << nlin << "\t"<< DDD<< "\t"<< rho1<< "\t"<<rr<< endl;


    NewtonPostProcess(u, f, info);
}

}
