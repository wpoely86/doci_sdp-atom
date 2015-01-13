#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <hdf5.h>

#include "LocalMinimizer.h"
#include "OptIndex.h"
#include "BoundaryPoint.h"
#include "PotentialReducation.h"

/**
 * You still need to set the max_angle, delta_angle, start_temp and
 * delta_temp after creating the object.
 * @param mol the molecular data to use
 */
simanneal::LocalMinimizer::LocalMinimizer(const CheMPS2::Hamiltonian &mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   OptIndex index(*ham);

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   energy = 0;
}

simanneal::LocalMinimizer::LocalMinimizer(CheMPS2::Hamiltonian &&mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   OptIndex index(*ham);

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   energy = 0;
}

simanneal::LocalMinimizer::~LocalMinimizer() = default;

/**
 * @return the real energy (calculated + nuclear repulsion)
 */
double simanneal::LocalMinimizer::get_energy() const
{
   return energy + ham->getEconst();
}

/**
 * Calculate the energy with the current
 * molecular data
 */
void simanneal::LocalMinimizer::calc_energy()
{
   method->BuildHam(*ham);
   method->Run();

   energy = method->getEnergy();
}

/**
 * Calculate the energy with the current
 * molecular data
 */
double simanneal::LocalMinimizer::calc_new_energy()
{
   orbtrans->fillHamCI(*ham);

   method->BuildHam(*ham);
   method->Run();

   return method->getEnergy();
}

simanneal::UnitaryMatrix& simanneal::LocalMinimizer::get_Optimal_Unitary()
{
   return orbtrans->get_unitary();
}

CheMPS2::Hamiltonian& simanneal::LocalMinimizer::getHam() const
{
   return *ham;
}

simanneal::OrbitalTransform& simanneal::LocalMinimizer::getOrbitalTf() const
{
   return *orbtrans;
}

doci2DM::Method& simanneal::LocalMinimizer::getMethod() const
{
   return *method;
}

/**
 * Return a reference to a PotentialReduction object. Only works if the
 * actual method is PotentialReduction. You should always wrap this in
 * a try/catch block and check for a std::bad_cast exception
 * @return PotentialReduction object of the current method
 */
doci2DM::PotentialReduction& simanneal::LocalMinimizer::getMethod_PR() const
{
   doci2DM::PotentialReduction &meth = dynamic_cast<doci2DM::PotentialReduction &> (*method);

   return meth;
}

/**
 * Return a reference to a BoundaryPoint object. Only works if the
 * actual method is BoundaryPoint. You should always wrap this in
 * a try/catch block and check for a std::bad_cast exception
 * @return BoundaryPoint object of the current method
 */
doci2DM::BoundaryPoint& simanneal::LocalMinimizer::getMethod_BP() const
{
   doci2DM::BoundaryPoint &meth = dynamic_cast<doci2DM::BoundaryPoint &> (*method);

   return meth;
}

void simanneal::LocalMinimizer::UseBoundaryPoint()
{
   method.reset(new doci2DM::BoundaryPoint(*ham));
}

void simanneal::LocalMinimizer::UsePotentialReduction()
{
   method.reset(new doci2DM::PotentialReduction(*ham));
}

std::vector< std::tuple<int,int,double,double> > simanneal::LocalMinimizer::scan_orbitals()
{
   const auto& ham2 = *ham;
   std::function<double(int,int)> getT = [&ham2] (int a, int b) -> double { return ham2.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&ham2]  (int a, int b, int c, int d) -> double { return ham2.getVmat(a,b,c,d); };

   std::vector< std::tuple<int,int,double,double> > pos_rotations;
   // worst case: c1 symmetry
   pos_rotations.reserve(ham->getL()*(ham->getL()-1)/2);

   for(int k_in=0;k_in<ham->getL();k_in++)
      for(int l_in=k_in+1;l_in<ham->getL();l_in++)
         if(ham->getOrbitalIrrep(k_in) == ham->getOrbitalIrrep(l_in))
         {
            auto found = method->getRDM().find_min_angle(k_in,l_in,0.3,getT,getV);

            if(!found.second)
               // we hit a maximum
               found = method->getRDM().find_min_angle(k_in,l_in,0.01,getT,getV);

            if(!found.second)
               // we're still stuck in a maximum, skip this!
               continue;

            double new_en = method->getRDM().calc_rotate(k_in,l_in,found.first,getT,getV);

            pos_rotations.push_back(std::make_tuple(k_in,l_in,found.first,new_en));
         }

   return pos_rotations;
}

/**
 * Do the local minimization
 */
void simanneal::LocalMinimizer::Minimize()
{
   bool converged = false;

   energy = method->getRDM().ddot(method->getHam());

   while(!converged)
   {
      auto list_rots = scan_orbitals();

      std::sort(list_rots.begin(), list_rots.end(),
            [](const std::tuple<int,int,double,double> & a, const std::tuple<int,int,double,double> & b) -> bool
            {
            return std::get<3>(a) < std::get<3>(b);
            });

      for(auto& elem: list_rots)
         std::cout << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t" << std::get<3>(elem) << std::endl;

      const auto& new_rot = list_rots[0];

      orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));

      auto new_energy = calc_new_energy();

      std::cout << "Rotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << std::get<3>(new_rot) << "  E = " << new_energy << std::endl;

      if(fabs(energy-new_energy)<1e-6)
         converged = true;

      energy = new_energy;
   }
}

/* vim: set ts=3 sw=3 expandtab :*/
