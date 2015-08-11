/* 
 * @BEGIN LICENSE
 *
 * Copyright (C) 2014-2015  Ward Poelmans
 *
 * This file is part of v2DM-DOCI.
 * 
 * v2DM-DOCI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * v2DM-DOCI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with v2DM-DOCI.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @END LICENSE
 */

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

      void WriteToFile(std::string filename) const;

      void ReadFromFile(std::string filename);

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
