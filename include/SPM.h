#ifndef SPM_H
#define SPM_H

#include <iostream>

#include "BlockStructure.h"

namespace doci2DM
{

class TPM;
class PHM;

class SPM: public BlockVector
{
   friend std::ostream &operator<<(std::ostream &output,SPM &spm);

   public:

      SPM(int L, int N);

      SPM(const TPM &);

      SPM(const SPM &) = default;

      SPM(SPM &&) = default;

      virtual ~SPM() = default;

      SPM& operator=(const SPM &) = default;

      SPM& operator=(SPM &&) = default;

      using BlockVector::operator=;

      using BlockVector::operator[];

      using BlockVector::operator();

      // do NOT override the operator()(int,int) from BlockVector
      double GetElement(int, int) const;

      int gN() const;

      int gL() const;

      int gn() const;

      void bar(double, const TPM &);

      void bar2(double, const TPM &);

      void bar3(double, const TPM &);

      void bar(double, const PHM &);

      void PrintSorted() const;

   private:

      //! number of particles
      int N;

      //! the size of the sp DOCI space (there are 2*L sp states)
      int L;
};

}

#endif /* SPM_H */

/* vim: set ts=3 sw=3 expandtab :*/
