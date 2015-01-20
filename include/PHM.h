#ifndef PHM_H
#define PHM_H

#include<memory>

#include "helpers.h"
#include "BlockStructure.h"

namespace doci2DM
{

class TPM;

class PHM: public BlockSymMatrix
{
   friend std::ostream &operator<<(std::ostream &output,PHM &phm);

   public:

      PHM(int, int);

      PHM(const PHM &) = default;

      PHM(PHM &&) = default;

      virtual ~PHM() = default;

      PHM& operator=(const PHM &) = default;

      PHM& operator=(PHM &&) = default;

      using BlockSymMatrix::operator=;

      using BlockSymMatrix::operator();

      double operator()(int a, int b, int c, int d) const;

      int gN() const;

      int gL() const;

      void G(const TPM &);

      SymMatrix Gimg(const TPM &) const;

      SymMatrix Gbuild() const;

      void sep_pm(BlockSymMatrix &, BlockSymMatrix &);

      void sqrt(int);

      void invert();

      void L_map(const BlockSymMatrix &, const BlockSymMatrix &);

   private:

      void constr_lists(int L);

      int L;

      int N;

      //! table translating single particles indices to two particle indices
      static std::unique_ptr<helpers::tmatrix<unsigned int>> s2ph;

      //! table translating two particles indices to single particle indices
      static std::unique_ptr<helpers::tmatrix<unsigned int>> ph2s;

      //! table translating single particles indices to the correct 2x2 block
      static std::unique_ptr<helpers::tmatrix<unsigned int>> s2b;

      //! table translating the block index to the single particle indices
      static std::unique_ptr<helpers::tmatrix<unsigned int>> b2s;
};

}

#endif /* PHM_H */

/*  vim: set ts=3 sw=3 expandtab :*/ 
