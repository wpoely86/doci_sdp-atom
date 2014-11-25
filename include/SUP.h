#ifndef SUP_H
#define SUP_H

#include <iostream>
#include <memory>
#include <string>

#include "include.h"

namespace doci2DM
{

class SUP
{
   friend std::ostream &operator<<(std::ostream &output,doci2DM::SUP &sup);

   public:

      SUP(int L, int N);

      SUP(const SUP &);

      SUP(SUP &&);

      virtual ~SUP() = default;

      SUP& operator=(const SUP &);

      SUP& operator=(SUP &&);

      SUP& operator=(double);

      SUP& operator+=(const SUP &);

      SUP& operator-=(const SUP &);

      SUP& operator*=(double);

      SUP& operator/=(double);

      void dscal(double);

      int gN() const;

      int gL() const;

      TPM const & getI() const;

      TPM & getI();

      TPM const & getQ() const;

      TPM & getQ();

      PHM const & getG() const;

      PHM & getG();

      void invert();

      void fill(const TPM &);

      void sqrt(int);

      void L_map(const SUP &, const SUP &);

      int gnr() const;

      double ddot(const SUP &) const;

      void daxpy(double, const SUP &);

      void sep_pm(SUP &, SUP &);

      void init_S(const Lineq &);

   private:
      //! number of particles
      int N;

      //! the size of the sp DOCI space (there are 2*L sp states)
      int L;

      //! the RDM matrix
      std::unique_ptr<TPM> I;

      //! the Q matrix
      std::unique_ptr<TPM> Q;

      //! the G matrix
      std::unique_ptr<PHM> G;
};

}

#endif /* SUP_H */

/* vim: set ts=3 sw=3 expandtab :*/
