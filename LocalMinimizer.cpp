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

   conv_crit = 1e-6;
   conv_steps = 100;
}

simanneal::LocalMinimizer::LocalMinimizer(CheMPS2::Hamiltonian &&mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   OptIndex index(*ham);

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   energy = 0;

   conv_crit = 1e-6;
   conv_steps = 100;
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

double simanneal::LocalMinimizer::calc_new_energy(const CheMPS2::Hamiltonian &new_ham)
{
   method->BuildHam(new_ham);
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
   auto start = std::chrono::high_resolution_clock::now();

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

            // skip angles larger than Pi/2
            if(fabs(found.first)>M_PI)
               continue;

            double new_en = method->getRDM().calc_rotate(k_in,l_in,found.first,getT,getV);

            assert(found.second && "Shit, maximum!");

            pos_rotations.push_back(std::make_tuple(k_in,l_in,found.first,new_en));
         }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Orbital scanning took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   assert(pos_rotations.size()>0);

   return pos_rotations;
}

/**
 * Do the local minimization
 */
void simanneal::LocalMinimizer::Minimize()
{
   int converged = 0;

   // first run
   method->Run();

   energy = method->getRDM().ddot(method->getHam());

   auto start = std::chrono::high_resolution_clock::now();

   std::pair<int,int> prev_pair(0,0);

   int iters = 0;

   while(converged<conv_steps)
   {
      auto list_rots = scan_orbitals();

      std::sort(list_rots.begin(), list_rots.end(),
            [](const std::tuple<int,int,double,double> & a, const std::tuple<int,int,double,double> & b) -> bool
            {
            return std::get<3>(a) < std::get<3>(b);
            });

      for(auto& elem: list_rots)
         std::cout << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t" << std::get<3>(elem)+ham->getEconst() << "\t" << std::get<2>(elem) << std::endl;

      std::pair<int,int> tmp = std::make_pair(std::get<0>(list_rots[0]), std::get<1>(list_rots[0]));
      int idx = 0;

      // don't do the same pair twice in a row
      if(tmp==prev_pair)
         idx++;

      const auto& new_rot = list_rots[idx];
      prev_pair = std::make_pair(std::get<0>(new_rot), std::get<1>(new_rot));

      assert(ham->getOrbitalIrrep(std::get<0>(new_rot)) == ham->getOrbitalIrrep(std::get<1>(new_rot)));
      // do Jacobi rotation twice: once for the Hamiltonian data and once for the Unitary Matrix
      orbtrans->DoJacobiRotation(*ham, std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));
      orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));

      double new_energy = calc_new_energy(*ham);

      if(fabs(energy-new_energy)<conv_crit)
         converged++;
      else 
         converged = 0;

      std::cout << iters << " (" << converged << ")\tRotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << std::get<3>(new_rot)+ham->getEconst() << "  E = " << new_energy+ham->getEconst() << "\t" << fabs(energy-new_energy) << std::endl;


      energy = new_energy;

      std::stringstream h5_name;
      h5_name << getenv("SAVE_H5_PATH") << "/unitary-" << iters << ".h5";
      orbtrans->get_unitary().saveU(h5_name.str());

      h5_name.str("");
      h5_name << getenv("SAVE_H5_PATH") << "/ham-" << iters << ".h5";
      ham->save2(h5_name.str());

      h5_name.str("");
      h5_name << getenv("SAVE_H5_PATH") << "/rdm-" << iters << ".h5";
      method->getRDM().WriteToFile(h5_name.str());

      doci2DM::BoundaryPoint *obj_bp = dynamic_cast<doci2DM::BoundaryPoint *> (method.get());
      if(obj_bp)
      {
         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/X-" << iters << ".h5";
         obj_bp->getX().WriteToFile(h5_name.str());

         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/Z-" << iters << ".h5";
         obj_bp->getZ().WriteToFile(h5_name.str());

         // if energy goes up instead of down or we have done 20 iters
         // then restart again from scratch
         if( ((std::get<3>(new_rot) - energy) < 0) || iters%20==0)
         {
            obj_bp->getX() = 0;
            obj_bp->getZ() = 0;
         }
      }

      iters++;
   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Minimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   std::stringstream h5_name;
   h5_name << getenv("SAVE_H5_PATH") << "/optimale-uni.h5";
   get_Optimal_Unitary().saveU(h5_name.str());
}

double simanneal::LocalMinimizer::get_conv_crit() const
{
   return conv_crit;
}

void simanneal::LocalMinimizer::set_conv_crit(double crit)
{
   conv_crit = crit;
}

void simanneal::LocalMinimizer::set_conv_steps(int steps)
{
   conv_steps = steps;
}

/* vim: set ts=3 sw=3 expandtab :*/
