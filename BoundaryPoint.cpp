#include <functional>
#include "BoundaryPoint.h"
#include "Hamiltonian.h"

using CheMPS2::Hamiltonian;
using doci2DM::BoundaryPoint;

BoundaryPoint::BoundaryPoint(const CheMPS2::Hamiltonian &hamin)
{
   N = Tools::getNumberOfParticles("CheMPS2_Ham_parent.h5");
   L = hamin.getL();
   nuclrep = hamin.getEconst();

   ham.reset(new TPM(L,N));

   X.reset(new SUP(L,N));
   Z.reset(new SUP(L,N));

   lineq.reset(new Lineq(L,N));

   BuildHam(hamin);

   // some default values
   sigma = 1.0;

   tol_PD = 3.0e-6;
   tol_en = 1.0e-3;

   mazzy = 1.0;

   max_iter = 5;
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
void BoundaryPoint::Run()
{
   TPM ham_copy(*ham);

   //only traceless hamiltonian needed in program.
   ham_copy.Proj_E(*lineq);

   //primal
   SUP X(L,N);

   //dual
   SUP Z(L,N);

   //Lagrange multiplier
   SUP V(L,N);

   //just dubya
   SUP W(L,N);

   SUP u_0(L,N);

   //little help
   TPM hulp(L,N);

   u_0.init_S(*lineq);

   X = 0.0;
   Z = 0.0;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   int iter_dual,iter_primal(0);

   int tot_iter = 0;

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
         SUP B(Z);

         B -= u_0;

         B.daxpy(1.0/sigma,X);

         TPM b(L,N);

         b.collaps(B, *lineq);

         b.daxpy(-1.0/sigma,ham_copy);

         hulp.InverseS(b, *lineq);

         hulp.Proj_E(*lineq);

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-mazzy/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v(L,N);

         v.collaps(V, *lineq);

         v -= ham_copy;

         D_conv = sqrt(v.ddot(v));
     }

      //update primal:
      X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z.getI().ddot(ham_copy) + X.ddot(u_0);

      std::cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z.getI().ddot(*ham) + nuclrep << "\t" << Z.getI().S_2() << std::endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;
   }

   energy = ham->ddot(Z.getI());

   std::cout << std::endl;
   std::cout << "Energy: " << ham->ddot(Z.getI()) + nuclrep << std::endl;
   std::cout << "Trace: " << Z.getI().trace() << std::endl;
   std::cout << "pd gap: " << Z.ddot(X) << std::endl;
   std::cout << "S^2: " << Z.getI().S_2() << std::endl;
   std::cout << "dual conv: " << D_conv << std::endl;
   std::cout << "primal conv: " << P_conv << std::endl;

   std::cout << std::endl;
   std::cout << "total nr of iterations = " << tot_iter << std::endl;
}

/**
 * @return the energy without the nuclear replusion part
 */
double BoundaryPoint::getEnergy() const
{
    return energy;
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

/* vim: set ts=3 sw=3 expandtab :*/