#include <functional>
#include "PotentialReducation.h"
#include "Hamiltonian.h"

using CheMPS2::Hamiltonian;
using doci2DM::PotentialReduction;

PotentialReduction::PotentialReduction(const CheMPS2::Hamiltonian &hamin)
{
   N = hamin.getNe();
   L = hamin.getL();
   nuclrep = hamin.getEconst();

   ham.reset(new TPM(L,N));

   rdm.reset(new TPM(L,N));

   lineq.reset(new Lineq(L,N));

   BuildHam(hamin);

   // some default values
   tolerance = 1.0e-5;
   target = 1e-12;
   reductionfac = 1.0/2.0;
}

PotentialReduction::PotentialReduction(const PotentialReduction &orig)
{
   N = orig.N;
   L = orig.L;
   nuclrep = orig.nuclrep;

   ham.reset(new TPM(*orig.ham));

   rdm.reset(new TPM(*orig.rdm));

   lineq.reset(new Lineq(*orig.lineq));

   // some default values
   tolerance = orig.tolerance;
   target = orig.target;
   reductionfac = orig.reductionfac;
   energy = orig.energy;
}

PotentialReduction& PotentialReduction::operator=(const PotentialReduction &orig)
{
   N = orig.N;
   L = orig.L;
   nuclrep = orig.nuclrep;

   (*ham) = *orig.ham;

   (*rdm) = *orig.rdm;

   (*lineq) = *orig.lineq;

   // some default values
   tolerance = orig.tolerance;
   target = orig.target;
   reductionfac = orig.reductionfac;
   energy = orig.energy;

   return *this;
}

PotentialReduction* PotentialReduction::Clone() const
{
   return new PotentialReduction(*this);
}

PotentialReduction* PotentialReduction::Move()
{
   return new PotentialReduction(std::move(*this));
}

/**
 * Build the new reduced hamiltonian based on the integrals
 * in ham
 * @param ham the integrals to use
 */
void PotentialReduction::BuildHam(const CheMPS2::Hamiltonian &hamin)
{
   std::function<double(int,int)> getT = [&hamin] (int a, int b) -> double { return hamin.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&hamin] (int a, int b, int c, int d) -> double { return hamin.getVmat(a,b,c,d); };

   ham->ham(getT, getV);

   norm_ham = std::sqrt(ham->ddot(*ham));
   (*ham) /= norm_ham;
}

/**
 * Do an actual calculation: calcalute the energy of the
 * reduced hamiltonian in ham
 */
void PotentialReduction::Run()
{
   rdm->init(*lineq);

   int tot_iter = 0.0;

   double t = 1.0;
   int iter = 0;

   TPM backup_rdm(*rdm);

   //outer iteration: scaling of the potential barrier
   while(t > target)
   {
      std::cout << iter << "\t" << t << "\t" << rdm->getMatrices().trace() << "\t" << rdm->getVectors().trace() << "\t" << rdm->ddot(*ham)*norm_ham + nuclrep << "\t" << rdm->S_2() << std::endl;

      double convergence = 1.0;
      iter++;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance)
      {
         tot_iter++;

         SUP P(L,N);

         P.fill(*rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(L,N);

         grad.constr_grad(t,P,*ham,*lineq);

         //dit wordt de stap:
         TPM delta(L,N);

         //los het hessiaan stelsel op:
         std::cout << delta.solve(t,P,grad,*lineq) << std::endl;

         //line search
         double a = delta.line_search(t,P,*ham);

         //rdm += a*delta;
         rdm->daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);
      }

      std::cout << std::endl;
      t *= reductionfac;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < target)
         tolerance = target;

      //extrapolatie:
      TPM extrapol(*rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = *rdm;

      double a = extrapol.line_search(t,*rdm,*ham);

      rdm->daxpy(a,extrapol);
   } 

   energy = norm_ham*ham->ddot(*rdm);

   std::cout << std::endl;
   std::cout << "Energy: " << getFullEnergy() << std::endl;
   std::cout << "Trace: " << rdm->trace() << std::endl;
   std::cout << "S^2: " << rdm->S_2() << std::endl;

   std::cout << std::endl;
   std::cout << "total nr of iterations = " << tot_iter << std::endl;
}

/**
 * @return the full energy (with the nuclear replusion part)
 */
double PotentialReduction::getFullEnergy() const
{
    return energy + nuclrep;
}

void PotentialReduction::set_target(double tar)
{
    this->target = tar;
}

void PotentialReduction::set_tolerance(double tol)
{
    this->tolerance = tol;
}

void PotentialReduction::set_reduction(double red)
{
    this->reductionfac = red;
}

doci2DM::TPM& PotentialReduction::getRDM() const
{
    return (*rdm);
}

doci2DM::Lineq& PotentialReduction::getLineq() const
{
    return (*lineq);
}

/* vim: set ts=3 sw=3 expandtab :*/
