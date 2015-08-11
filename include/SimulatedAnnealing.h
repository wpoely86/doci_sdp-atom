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

#ifndef SIM_ANNEAL_H
#define SIM_ANNEAL_H

#include <random>

#include "Method.h"
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"

namespace doci2DM {
class PotentialReduction;
class BoundaryPoint;
}

namespace simanneal
{ 

class SimulatedAnnealing
{
   public:
      SimulatedAnnealing(const CheMPS2::Hamiltonian &);

      SimulatedAnnealing(CheMPS2::Hamiltonian &&);

      virtual ~SimulatedAnnealing();

      bool accept_function(double);

      void optimize();

      void optimize_mpi();

      double calc_new_energy();

      void calc_energy();

      double get_energy() const;

      void Set_max_angle(double);

      void Set_delta_angle(double);

      void Set_start_temp(double);

      void Set_delta_temp(double);

      simanneal::UnitaryMatrix& get_Optimal_Unitary();

      CheMPS2::Hamiltonian& getHam() const;

      OrbitalTransform& getOrbitalTf() const;

      doci2DM::Method& getMethod() const;

      void UseBoundaryPoint();

      void UsePotentialReduction();

      doci2DM::PotentialReduction& getMethod_PR() const;

      doci2DM::BoundaryPoint& getMethod_BP() const;

   private:

      std::unique_ptr<doci2DM::Method> method;

      //! Holds the current hamiltonian
      std::unique_ptr<CheMPS2::Hamiltonian> ham;

      //! the actual orbital transform
      std::unique_ptr<simanneal::OrbitalTransform> orbtrans;

      //! the current unitary
      std::unique_ptr<simanneal::UnitaryMatrix> opt_unitary;

      //! energy of current iteration
      double energy;
      //! the start temperatur
      double start_temp;
      //! the change in temperatur between steps
      double delta_temp;
      //! the change in the angle allows in steps
      double delta_angle;
      //! the maximum allows angle
      double max_angle;
      //! the number of steps done
      unsigned int steps;
      //! max number of steps
      unsigned int max_steps;

      //! current temperature
      double cur_temp;

      //! number of unaccepted steps
      unsigned int unaccepted;

      //! bool to indicate if we should quite the optimalisation
      bool stop_running;

      //! the real random input (hopefully)
      std::random_device rd;
      //! our pseudo-random generator
      std::mt19937_64 mt;
};

}

#endif /* SIM_ANNEAL_H */

/* vim: set ts=3 sw=3 expandtab :*/
