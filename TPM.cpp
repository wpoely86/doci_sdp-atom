#include <assert.h>

#include "TPM.h"

/**
 * Create a object in the Two Particle space
 * @param L the number of levels (the sp space has size of 2*L)
 * @param N the number of particles
 */
TPM::TPM(int L, int N): Container(1,1)
{
   this->N = N;
   this->L = L;
   n = L*(2*L-1);

   // the LxL block with degen = 1
   setMatrixDim(0, L, 1);
   // the L/2*(L-1) vector with degen = 4
   setVectorDim(0, (L*(L-1))/2, 4);
}

//TPM::TPM(const TPM &orig): Container(orig)
//{
//   N = orig.N;
//   L = orig.L;
//   n = orig.n;
//}
//
//TPM::TPM(TPM &&orig): Container(orig)
//{
//   N = orig.N;
//   L = orig.L;
//   n = orig.n;
//}

TPM::~TPM()
{
}

int TPM::gN() const
{
   return N;
}

int TPM::gL() const
{
   return L;
}

int TPM::gn() const
{
   return n;
}

std::ostream &operator<<(std::ostream &output,TPM &tpm)
{
   output << "Block: " << std::endl;
   output << tpm.getMatrix(0) << std::endl;
   output << std::endl;
   output << "Vector (4x): " << std::endl;
   output << tpm.getVector(0) << std::endl;

   return output;
}

/*  vim: set ts=3 sw=3 expandtab :*/
