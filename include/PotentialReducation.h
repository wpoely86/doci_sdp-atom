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

#ifndef POTENTIAL_REDUCTION_H
#define POTENTIAL_REDUCTION_H

#include "include.h"

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{

class PotentialReduction: public Method
{
   public:

      PotentialReduction(const CheMPS2::Hamiltonian &);

      PotentialReduction(const TPM &);

      PotentialReduction(const PotentialReduction &);

      PotentialReduction(PotentialReduction &&) = default;

      virtual ~PotentialReduction() = default;

      PotentialReduction* Clone() const;

      PotentialReduction* Move();

      PotentialReduction& operator=(PotentialReduction &&) = default;

      PotentialReduction& operator=(const PotentialReduction &);

      void BuildHam(const CheMPS2::Hamiltonian &);

      void BuildHam(const TPM &);

      unsigned int Run();

      double getFullEnergy() const;

      void set_target(double);

      void set_tolerance(double);

      void set_reduction(double);

      TPM& getRDM() const;

      TPM& getHam() const;

      Lineq& getLineq() const;

      double evalEnergy() const;

      bool FullyConverged() const;

   private:

      std::unique_ptr<TPM> ham;

      std::unique_ptr<TPM> rdm;

      std::unique_ptr<Lineq> lineq;

      double nuclrep;

      double norm_ham;

      double target;

      double tolerance;

      double reductionfac;

      double t;
};

}

#endif /* POTENTIAL_REDUCTION_H */

/* vim: set ts=3 sw=3 expandtab :*/
