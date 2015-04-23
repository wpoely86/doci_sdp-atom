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
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @END LICENSE
 */

#ifndef SPM_H
#define SPM_H

#include <iostream>

#include "BlockStructure.h"

namespace simanneal { class OptIndex; }

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

      std::vector<double> Particlesperirrep(const simanneal::OptIndex &index, bool print=true) const;

   private:

      //! number of particles
      int N;

      //! the size of the sp DOCI space (there are 2*L sp states)
      int L;
};

}

#endif /* SPM_H */

/* vim: set ts=3 sw=3 expandtab :*/
