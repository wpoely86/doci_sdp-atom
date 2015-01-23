#ifndef LOCALMINIMIZER_H
#define LOCALMINIMIZER_H

#include <memory>
#include <vector>
#include <tuple>

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

      void Minimize();

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

      std::vector< std::tuple<int,int,double,double> > scan_orbitals();

      double get_conv_crit() const;

      void set_conv_crit(double);

      void set_conv_steps(int);

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
};

}

#endif /* LOCALMINIMIZER_H */

/*  vim: set ts=3 sw=3 expandtab :*/
