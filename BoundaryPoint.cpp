#include <fstream>
#include <chrono>
#include <functional>
#include "BoundaryPoint.h"
#include "Hamiltonian.h"

using CheMPS2::Hamiltonian;
using doci2DM::BoundaryPoint;

BoundaryPoint::BoundaryPoint(const CheMPS2::Hamiltonian &hamin)
{
   N = hamin.getNe();
   L = hamin.getL();
   nuclrep = hamin.getEconst();

   ham.reset(new TPM(L,N));

   X.reset(new SUP(L,N));
   Z.reset(new SUP(L,N));

   useprevresult = false;
   (*X) = 0.0;
   (*Z) = 0.0;

   lineq.reset(new Lineq(L,N));

   BuildHam(hamin);

   // some default values
   sigma = 1.0;

   tol_PD = 3.0e-6;
   tol_en = 1.0e-3;

   mazzy = 1.0;

   max_iter = 5;

   avg_iters = 100000000; // first step we don't really limited anything
   iters = 0;
   runs = 0;
}

BoundaryPoint::BoundaryPoint(const BoundaryPoint &orig)
{
   N = orig.N;
   L = orig.L;
   nuclrep = orig.nuclrep;

   ham.reset(new TPM(*orig.ham));

   X.reset(new SUP(*orig.X));
   Z.reset(new SUP(*orig.Z));

   lineq.reset(new Lineq(*orig.lineq));

   useprevresult = orig.useprevresult;

   sigma = orig.sigma;;

   tol_PD = orig.tol_PD;
   tol_en = orig.tol_en;

   mazzy = orig.mazzy;

   max_iter = orig.max_iter;

   energy = orig.energy;

   avg_iters = orig.avg_iters;
   iters = orig.iters;
   runs = orig.runs;
}

BoundaryPoint& BoundaryPoint::operator=(const BoundaryPoint &orig)
{
   N = orig.N;
   L = orig.L;
   nuclrep = orig.nuclrep;

   (*ham) = *orig.ham;

   (*X) = *orig.X;
   (*Z) = *orig.Z;

   (*lineq) = *orig.lineq;

   useprevresult = orig.useprevresult;

   sigma = orig.sigma;;

   tol_PD = orig.tol_PD;
   tol_en = orig.tol_en;

   mazzy = orig.mazzy;

   max_iter = orig.max_iter;

   energy = orig.energy;

   avg_iters = orig.avg_iters;
   iters = orig.iters;
   runs = orig.runs;

   return *this;
}

BoundaryPoint* BoundaryPoint::Clone() const
{
   return new BoundaryPoint(*this);
}

BoundaryPoint* BoundaryPoint::Move()
{
   return new BoundaryPoint(std::move(*this));
}

/**
 * Build the new reduced hamiltonian based on the integrals
 * in ham
 * @param ham the integrals to use
 */
void BoundaryPoint::BuildHam(const CheMPS2::Hamiltonian &hamin)
{
   std::function<double(int,int)> getT = [&hamin] (int a, int b) -> double { return hamin.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&hamin] (int a, int b, int c, int d) -> double { return hamin.getVmat(a,b,c,d); };

   ham->ham(getT, getV);
}

/**
 * Do an actual calculation: calcalute the energy of the
 * reduced hamiltonian in ham
 */
unsigned int BoundaryPoint::Run()
{
   TPM ham_copy(*ham);

   //only traceless hamiltonian needed in program.
   ham_copy.Proj_E(*lineq);

   if(!useprevresult)
   {
      (*X) = 0;
      (*Z) = 0;
   }

   //Lagrange multiplier
   SUP V(L,N);

   //just dubya
   SUP W(L,N);

   SUP u_0(L,N);

   //little help
   TPM hulp(L,N);

   u_0.init_S(*lineq);

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   unsigned int iter_dual(0),iter_primal(0);

   unsigned int tot_iter = 0;

   std::ostream* fp = &std::cout;
   std::ofstream fout;
   if(!outfile.empty())
   {
      fout.open(outfile, std::ios::out | std::ios::app);
      fp = &fout;
   }
   std::ostream &out = *fp;
   out.precision(10);
   out.unsetf(std::ios_base::floatfield);

   auto start = std::chrono::high_resolution_clock::now();

   while(P_conv > tol_PD || D_conv > tol_PD) // || fabs(convergence) > tol_en)
   {
      ++iter_primal;

      D_conv = 1.0;

      iter_dual = 0;

      while(D_conv > tol_PD  && iter_dual <= max_iter)
      {
         ++tot_iter;

         ++iter_dual;

         //solve system
         SUP B(*Z);

         B -= u_0;

         B.daxpy(1.0/sigma,*X);

         TPM b(L,N);

         b.collaps(B, *lineq);

         b.daxpy(-1.0/sigma,ham_copy);

         hulp.InverseS(b, *lineq);

         hulp.Proj_E(*lineq);

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-mazzy/sigma,*X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(*Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v(L,N);

         v.collaps(V, *lineq);

         v -= ham_copy;

         D_conv = sqrt(v.ddot(v));
     }

      //update primal:
      *X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= *Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z->getI().ddot(ham_copy) + X->ddot(u_0);

      if(do_output && iter_primal%500 == 0)
         out << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z->getI().ddot(*ham) + nuclrep << "\t" << Z->getI().S_2() << std::endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

      if(iter_primal>avg_iters*5)
         break;
   }

   auto end = std::chrono::high_resolution_clock::now();

   if(iter_primal>avg_iters*5)
      energy = 1e24; // something big so we're sure the step will be rejected
   else
   {
      energy = ham->ddot(Z->getI());
      runs++;
      iters += iter_primal;
      avg_iters = iters/runs;
   }

   out << std::endl;
   out << "Energy: " << ham->ddot(Z->getI()) + nuclrep << std::endl;
   out << "Trace: " << Z->getI().trace() << std::endl;
   out << "pd gap: " << Z->ddot(*X) << std::endl;
   out << "S^2: " << Z->getI().S_2() << std::endl;
   out << "dual conv: " << D_conv << std::endl;
   out << "primal conv: " << P_conv << std::endl;
   out << "Runtime: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;
   out << "Primal iters: " << iter_primal << std::endl;
   out << "avg primal iters: " << avg_iters << std::endl;

   out << std::endl;
   out << "total nr of iterations = " << tot_iter << std::endl;

   if(!outfile.empty())
      fout.close();

   return iter_primal;
}

/**
 * @return the full energy (with the nuclear replusion part)
 */
double BoundaryPoint::getFullEnergy() const
{
    return energy + nuclrep;
}

void BoundaryPoint::set_tol_PD(double tol)
{
    this->tol_PD = tol;
}

void BoundaryPoint::set_tol_en(double tol)
{
    this->tol_en = tol;
}

void BoundaryPoint::set_mazzy(double maz)
{
    this->mazzy = maz;
}

void BoundaryPoint::set_sigma(double sig)
{
    this->sigma = sig;
}

void BoundaryPoint::set_max_iter(int iters)
{
    this->max_iter = iters;
}

doci2DM::SUP& BoundaryPoint::getX() const
{
    return (*X);
}

doci2DM::SUP& BoundaryPoint::getZ() const
{
    return (*Z);
}

doci2DM::Lineq& BoundaryPoint::getLineq() const
{
    return (*lineq);
}


doci2DM::TPM& BoundaryPoint::getRDM() const
{
   return Z->getI();
}


void BoundaryPoint::set_use_prev_result(bool new_val)
{
   useprevresult = new_val;
}

/* vim: set ts=3 sw=3 expandtab :*/
