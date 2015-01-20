#include "include.h"
#include "lapack.h"

using namespace doci2DM;

/**
 * Create SPM matrix (vector)
 * @param L the number of orbitals (spinorbitals/2)
 * @param N the number of particles
 */
SPM::SPM(int L, int N): BlockVector(1)
{
   this->L = L;
   this->N = N;

   setDim(0,L,2);
}

/**
 * Construct the 1DM from the 2DM in the most 
 * general way
 * @param tpm the 2DM to use
 */
SPM::SPM(const TPM &tpm): BlockVector(1)
{
   this->L = tpm.gL();
   this->N = tpm.gN();

   setDim(0,L,2);

   bar3(1.0/(N-1.0), tpm);
}

int SPM::gN() const
{
   return N;
}

int SPM::gL() const
{
   return L;
}

int SPM::gn() const
{
   return 2*L;
}

/**
 * Get element a,b from the SPM
 * This does not override the operator()(int,int)
 * from BlockVector
 * @param a the first sp index
 * @param b the second sp index
 * @return the requested element
 */
double SPM::GetElement(int a, int b) const
{
   if(a!=b)
      return 0;

   return (*this)(0,a%L);
}

/**
 * contruct the 1DM from the block in the 2DM
 * @param scal the scaling factor
 * @param tpm the 2DM to use
 */
void SPM::bar(double scal, const TPM &tpm)
{
   for(int i=0;i<L;i++)
      (*this)(0,i) = scal * tpm(0,i,i);

//   int incx = L+1;
//   int incy = 1;
//   dcopy_(&L, tpm.getMatrix(0).gMatrix(), &incx, (*this)[0].gVector(), &incy);
}

/**
 * Construct the 1DM using the vector part of the 2DM
 * @param scal the scaling factor
 * @param tpm the 2DM to use
 */
void SPM::bar2(double scal, const TPM &tpm)
{
   for(int i=0;i<L;i++)
   {
      (*this)(0,i) = 0;

      for(int a=0;a<L;a++)
         (*this)(0,i) += tpm.getDiag(i,a);

      (*this)(0,i) *= scal;
   }
}

/**
 * Calculate the 1DM in the most general way for DOCI.
 * This must be used in the real Q image
 * @param scal the scaling factor
 * @param tpm the 2DM to use
 */
void SPM::bar3(double scal, const TPM &tpm)
{
   for(int i=0;i<L;i++)
   {
      (*this)(0,i) = tpm.getMatrix(0)(i,i);

      for(int a=0;a<L;a++)
         (*this)(0,i) += 2 * tpm.getDiag(i,a);

      (*this)(0,i) *= scal;
   }
}


namespace doci2DM {
std::ostream &operator<<(std::ostream &output,SPM &spm)
{
   output << "SPM (2x): " << std::endl;

   for(int i=0;i<spm.gL();i++)
      output << i << "\t|\t" << i << "  " << i << "\t\t" << spm(0,i) << std::endl;

   return output;
}
}

void SPM::bar(double scal, const PHM &phm)
{
   for(int i=0;i<L;i++)
   {
      (*this)(0,i) = 0;

      for(int a=0;a<L;a++)
         (*this)(0,i) += phm(i,a,i,a)+phm(i,a+L,i,a+L);

      (*this)(0,i) *= scal;
   }
}

/*  vim: set ts=3 sw=3 expandtab :*/
