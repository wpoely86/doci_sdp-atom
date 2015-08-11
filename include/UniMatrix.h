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

#ifndef UNITARY_MATRIX_H
#define UNITARY_MATRIX_H

#include <iostream>

#include "BlockStructure.h"

namespace simanneal { class OptIndex; }

namespace doci2DM
{

/**
 * This class is for real unitary (block) matrices.
 */
class UnitaryMatrix: public BlockMatrix
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de Matrix you want to print
    */
   friend std::ostream &operator<<(std::ostream &,const doci2DM::UnitaryMatrix &);

   public:

      UnitaryMatrix(const simanneal::OptIndex &);

      UnitaryMatrix(int);

      UnitaryMatrix(const UnitaryMatrix &) = default;

      UnitaryMatrix(UnitaryMatrix &&) = default;

      virtual ~UnitaryMatrix() = default;

      UnitaryMatrix &operator=(const UnitaryMatrix &) = default;

      UnitaryMatrix &operator=(UnitaryMatrix &&) = default;

      using BlockMatrix::operator=;

   private:

      std::unique_ptr<simanneal::OptIndex> index;

};

}

#endif /* UNITARY_MATRIX_H */

/* vim: set ts=3 sw=3 expandtab :*/
