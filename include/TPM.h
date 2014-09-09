#ifndef TPM_H
#define TPM_H

#include <iostream>

#include "Container.h"

class TPM: public Container
{
   friend std::ostream &operator<<(std::ostream &output,TPM &tpm);

   public:

      TPM(int L, int N);

      TPM(const TPM &) = default;

      TPM(TPM &&) = default;

      virtual ~TPM();

      TPM& operator=(const TPM &) = default;

      TPM& operator=(TPM &&) = default;

      using Container::operator=;

      using Container::operator();

      int gN() const;

      int gL() const;

      int gn() const;

   private:

      //! number of particles
      int N;

      //! the size of the sp DOCI space (there are 2*L sp states)
      int L;

      //! dimension of the full TPM
      int n;
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
