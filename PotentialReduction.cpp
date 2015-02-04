#include <fstream>
#include <iomanip>
#include <chrono>
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
   target = 1e-11;
   reductionfac = 1.0/1.01;
   t = 1;
}

PotentialReduction::PotentialReduction(const TPM &hamin)
{
   N = hamin.gN();
   L = hamin.gL();
   nuclrep = 0;

   ham.reset(new TPM(hamin));

   rdm.reset(new TPM(L,N));

   lineq.reset(new Lineq(L,N));

   BuildHam(hamin);

   // some default values
   tolerance = 1.0e-5;
   target = 1e-12;
   reductionfac = 1.0/1.1;
   t = 1;
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
   t = orig.t;
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
   t = orig.t;

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
 * Copy the reduced hamiltonian from a TPM object
 * @param ham the TPM object to use
 */
void PotentialReduction::BuildHam(const TPM &hamin)
{
   (*ham) = hamin;

   norm_ham = std::sqrt(ham->ddot(*ham));
   (*ham) /= norm_ham;
}

/**
 * Do an actual calculation: calcalute the energy of the
 * reduced hamiltonian in ham
 */
unsigned int PotentialReduction::Run()
{
   rdm->init(*lineq);

   unsigned int tot_iter = 0;

   t = 1.0;
   int iter = 0;

   TPM backup_rdm(*rdm);

   std::ostream* fp = &std::cout;
   std::ofstream fout;
   if(!outfile.empty())
   {
      fout.open(outfile, std::ios::out | std::ios::app);
      fp = &fout;
    }
   std::ostream &out = *fp;
   out.precision(10);
   out.setf(std::ios::scientific | std::ios::fixed, std::ios_base::floatfield);

   auto start = std::chrono::high_resolution_clock::now();

   //outer iteration: scaling of the potential barrier
   while(t > target)
   {
      if(do_output)
         out << iter << "\t" << std::setw(16) << t << "\t" << std::setw(16) << rdm->getMatrices().trace() << "\t" << std::setw(16) << rdm->getVectors().trace() << "\t" << std::setw(16) << rdm->ddot(*ham)*norm_ham + nuclrep << "\t" << std::setw(16) << rdm->S_2() << std::endl;

      double convergence = 1.0;
      int cg_iters = 0;
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
         cg_iters = delta.solve(t,P,grad,*lineq);

         //line search
         double a = delta.line_search(t,P,*ham);

         //rdm += a*delta;
         rdm->daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

         if(do_output)
            out << cg_iters << "\t" << convergence << std::endl;

         if(cg_iters == -1)
            break;
      }

      if(do_output)
         out << std::endl;

      if(cg_iters == -1)
      {
         *rdm = backup_rdm;
         break;
      }

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

   auto end = std::chrono::high_resolution_clock::now();

   energy = norm_ham*ham->ddot(*rdm);

   out << std::endl;
   out << "Energy: " << getFullEnergy() << std::endl;
   out << "Trace: " << rdm->trace() << std::endl;
   out << "pd gap: " << t*rdm->gn() << std::endl;
   out << "S^2: " << rdm->S_2() << std::endl;
   out << "Runtime: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   out << std::endl;
   out << "total nr of iterations = " << tot_iter << std::endl;

   if(!outfile.empty())
      fout.close();

   return tot_iter;
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

doci2DM::TPM& PotentialReduction::getHam() const
{
   return *ham;
}

/**
 * Give the energy with the current rdm and ham
 * @param the newly evaluated energy
 */
double PotentialReduction::evalEnergy() const
{
   return norm_ham*ham->ddot(*rdm) + nuclrep;
}

/**
 * Check if last calculation was fully convergenced
 * @return true if all convergence critera are met, false otherwise
 */
bool PotentialReduction::FullyConverged() const
{
   return !(t > target);
}

/* vim: set ts=3 sw=3 expandtab :*/
