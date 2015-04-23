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

#ifndef LOCALMINIMIZER_H
#define LOCALMINIMIZER_H

#include <memory>
#include <vector>
#include <tuple>
#include <random>

#include "Method.h"
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"

namespace doci2DM {
   class PotentialReduction;
   class BoundaryPoint;
}

namespace simanneal {

class LocalMinimizer
{
   public:
      LocalMinimizer(const CheMPS2::Hamiltonian &);

      LocalMinimizer(CheMPS2::Hamiltonian &&);

      virtual ~LocalMinimizer();

      int Minimize(bool dist_choice=false, int start_iters=0);

      double get_energy() const;

      double calc_new_energy();

      double calc_new_energy(const CheMPS2::Hamiltonian &);

      void calc_energy();

      simanneal::UnitaryMatrix& get_Optimal_Unitary();

      CheMPS2::Hamiltonian& getHam() const;

      OrbitalTransform& getOrbitalTf() const;

      doci2DM::Method& getMethod() const;

      void UseBoundaryPoint();

      void UsePotentialReduction();

      doci2DM::PotentialReduction& getMethod_PR() const;

      doci2DM::BoundaryPoint& getMethod_BP() const;

      std::vector<std::tuple<int,int,double,double>> scan_orbitals();

      double get_conv_crit() const;

      void set_conv_crit(double);

      void set_conv_steps(int);

      int choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &);

      int Minimize_noOpt(double stopcrit);

      int Minimize_hybrid();

   private:

      //! criteria for convergence of the minimizer
      double conv_crit;

      double energy;

      //! number of steps in convergence area
      int conv_steps;

      std::unique_ptr<doci2DM::Method> method;

      //! Holds the current hamiltonian
      std::unique_ptr<CheMPS2::Hamiltonian> ham;

      //! the actual orbital transform
      std::unique_ptr<simanneal::OrbitalTransform> orbtrans;

      //! the current unitary
      std::unique_ptr<simanneal::UnitaryMatrix> opt_unitary;

      //! our pseudo-random generator
      std::mt19937 mt;

      //! only rotation within these irreps (if not empty)
      std::vector<int> allow_irreps;
};

}

#endif /* LOCALMINIMIZER_H */

/*  vim: set ts=3 sw=3 expandtab :*/
