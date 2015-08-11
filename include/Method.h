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

#ifndef METHOD_H
#define METHOD_H

#include <memory>

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{
class TPM;

class Method
{
   public:

      Method() { do_output = true; }

      virtual ~Method() = default;

      virtual unsigned int Run() = 0;

      virtual void BuildHam(const CheMPS2::Hamiltonian &) = 0;

      virtual void BuildHam(const TPM &) = 0;

      virtual Method* Clone() const = 0;

      virtual Method* Move() = 0;

      virtual TPM& getRDM() const = 0;

      virtual TPM& getHam() const = 0;

      virtual double evalEnergy() const = 0;

      double getEnergy() const { return energy; }

      virtual void set_output(bool out) { do_output = out; } 

      virtual void set_outfile(std::string filename) { outfile = filename; }

      virtual bool FullyConverged() const = 0;

   protected:

      int L;

      int N;

      double energy;

      bool do_output;

      std::string outfile;
};

}

#endif /* METHOD_H */

/* vim: set ts=3 sw=3 expandtab :*/
