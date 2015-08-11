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

#ifndef TOOLS_H
#define TOOLS_H

#include <string>

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{
   class TPM;

class Tools
{
   public:
      static int getspDimension(std::string filename);

      static int getNumberOfParticles(std::string filename);

      static double getNuclearRepulEnergy(std::string filename);

      static void scan_all(const TPM &rdm, const CheMPS2::Hamiltonian &ham);

      static void scan_all_bp(const TPM &rdm, const CheMPS2::Hamiltonian &ham);
};

}

#endif /* TOOLS_H */

/* vim: set ts=3 sw=3 expandtab :*/
