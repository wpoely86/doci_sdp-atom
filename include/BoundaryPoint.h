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

#ifndef BOUNDARY_POINT_H
#define BOUNDARY_POINT_H

#include "include.h"

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{

class BoundaryPoint: public Method
{
   public:

      BoundaryPoint(const CheMPS2::Hamiltonian &);

      BoundaryPoint(const TPM &);

      BoundaryPoint(const BoundaryPoint &);

      BoundaryPoint(BoundaryPoint &&) = default;

      virtual ~BoundaryPoint() = default;

      BoundaryPoint& operator=(const BoundaryPoint &);

      BoundaryPoint& operator=(BoundaryPoint &&) = default;

      BoundaryPoint* Clone() const;

      BoundaryPoint* Move();

      void BuildHam(const CheMPS2::Hamiltonian &);

      void BuildHam(const TPM &);

      unsigned int Run();

      double getFullEnergy() const;

      void set_tol_PD(double);

      void set_tol_en(double);

      void set_mazzy(double);

      void set_sigma(double);

      void set_max_iter(unsigned int);

      double get_tol_PD() const;

      SUP& getX() const;

      SUP& getZ() const;

      Lineq& getLineq() const;

      TPM& getRDM() const;

      TPM& getHam() const;

      void set_use_prev_result(bool);

      double evalEnergy() const;

      void ReturnHighWhenBailingOut(bool);

      void Reset_avg_iters();

      double get_P_conv() const;

      double get_D_conv() const;

      double get_convergence() const;

      bool FullyConverged() const;

      std::vector<double> energyperirrep(const CheMPS2::Hamiltonian &, bool print=false);

   private:

      std::unique_ptr<TPM> ham;

      std::unique_ptr<SUP> X;

      std::unique_ptr<SUP> Z;

      std::unique_ptr<Lineq> lineq;

      double nuclrep;

      double tol_PD, tol_en;

      double mazzy;

      double sigma;

      unsigned int max_iter;

      unsigned int avg_iters;

      unsigned int iters;
      unsigned int runs;

      bool useprevresult;

      //! when true, return very high value for the energy if the calculation takes too many iterations
      bool returnhigh;

      //! the 3 convergence criteria
      double D_conv, P_conv, convergence;
};

}

#endif /* BOUNDARY_POINT_H */

/* vim: set ts=3 sw=3 expandtab :*/
