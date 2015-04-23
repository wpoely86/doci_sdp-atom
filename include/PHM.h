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

#ifndef PHM_H
#define PHM_H

#include<memory>
#include<hdf5.h>

#include "helpers.h"
#include "BlockStructure.h"

namespace doci2DM
{

class TPM;

class PHM: public BlockMatrix
{
   friend std::ostream &operator<<(std::ostream &output,PHM &phm);

   public:

      PHM(int, int);

      PHM(const PHM &) = default;

      PHM(PHM &&) = default;

      virtual ~PHM() = default;

      PHM& operator=(const PHM &) = default;

      PHM& operator=(PHM &&) = default;

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int a, int b, int c, int d) const;

      const Matrix& getBlock(int a, int b) const;

      Matrix& getBlock(int a, int b);

      int gN() const;

      int gL() const;

      void G(const TPM &);

      Matrix Gimg(const TPM &) const;

      Matrix Gbuild() const;

      void sep_pm(BlockMatrix &, BlockMatrix &);

      void sqrt(int);

      void invert();

      void L_map(const BlockMatrix &, const BlockMatrix &);

      void WriteToFile(hid_t &group_id) const;

      void ReadFromFile(hid_t &group_id);

      void WriteFullToFile(hid_t &group_id) const;

   private:

      void constr_lists(int L);

      int L;

      int N;

      //! table translating single particles indices to two particle indices
      static std::unique_ptr<helpers::tmatrix<int>> s2ph;

      //! table translating two particles indices to single particle indices
      static std::unique_ptr<helpers::tmatrix<int>> ph2s;

      //! table translating single particles indices to the correct 2x2 block
      static std::unique_ptr<helpers::tmatrix<int>> s2b;

      //! table translating the block index to the single particle indices
      static std::unique_ptr<helpers::tmatrix<int>> b2s;
};

}

#endif /* PHM_H */

/*  vim: set ts=3 sw=3 expandtab :*/ 
