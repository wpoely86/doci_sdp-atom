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

#ifndef EIG_H
#define EIG_H

#include <iostream>

#include "include.h"

namespace doci2DM
{

class EIG: public BlockVector
{
   public:

      EIG(BlockMatrix &);

      EIG(Container &);

      EIG(SUP &);

      virtual ~EIG() = default;

      using BlockVector::operator=;

      using BlockVector::operator();

      using BlockVector::operator[];

      double min() const;

      double max() const;

      double lsfunc(double) const;
};

}

#endif /* EIG_H */

/* vim: set ts=3 sw=3 expandtab :*/
