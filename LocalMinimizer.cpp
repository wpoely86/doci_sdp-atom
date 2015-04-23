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

#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <hdf5.h>
#include <signal.h>
#include <cstring>

#include "LocalMinimizer.h"
#include "OptIndex.h"
#include "BoundaryPoint.h"
#include "PotentialReducation.h"

// if set, the signal has been given to stop the minimalisation
extern sig_atomic_t stopping_min;

/**
 * @param mol the molecular data to use
 */
simanneal::LocalMinimizer::LocalMinimizer(const CheMPS2::Hamiltonian &mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   energy = 0;

   conv_crit = 1e-6;
   conv_steps = 50;

   std::random_device rd;
   mt = std::mt19937(rd());

   // expect to be comma seperated list of allowed irreps
   char *irreps_env = getenv("v2DM_DOCI_ALLOWED_IRREPS");
   if(irreps_env && strlen(irreps_env) > 0)
   {
      std::string irreps_string = irreps_env;
      const std::string delim = ",";
      CheMPS2::Irreps syminfo(ham->getNGroup());

      auto start = 0U;
      auto end = irreps_string.find(delim);
      try
      { 
         while (true)
         {
            auto elem = irreps_string.substr(start, end - start);
            if(elem.empty())
               break;

            int cur_irrep = std::stoi(elem);

            if(cur_irrep >= 0 && cur_irrep < syminfo.getNumberOfIrreps())
               allow_irreps.push_back(cur_irrep);

            start = end + delim.length();

            if(end >= std::string::npos)
               break;

            end = irreps_string.find(delim, start);
         }
      } catch (std::exception& e) {
         std::cout << "Invalid value in v2DM_DOCI_ALLOWED_IRREPS" << std::endl;
      }

      std::sort(allow_irreps.begin(), allow_irreps.end());
      std::cout << "Allowed irreps: ";
      for(auto &elem: allow_irreps)
         std::cout << elem << " ";
      std::cout << std::endl;
   }
}

simanneal::LocalMinimizer::LocalMinimizer(CheMPS2::Hamiltonian &&mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   energy = 0;

   conv_crit = 1e-6;
   conv_steps = 50;

   std::random_device rd;
   mt = std::mt19937(rd());
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
            if(!allow_irreps.empty() && std::find(allow_irreps.begin(), allow_irreps.end(), ham->getOrbitalIrrep(k_in)) == allow_irreps.end() )
               continue;

            auto found = method->getRDM().find_min_angle(k_in,l_in,0.3,getT,getV);

            if(!found.second)
               // we hit a maximum
               found = method->getRDM().find_min_angle(k_in,l_in,0.01,getT,getV);

            if(!found.second)
               // we're still stuck in a maximum, skip this!
               continue;

            // skip angles larger than Pi/2
            if(fabs(found.first)>M_PI/2.0)
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
 * @param dist_choice if set to true, we use choose_orbitals to choose
 * which pair of orbitals to use (instead of the lowest one)
 * @param start_iters start number the iterations from this number (defaults to 0)
 */
int simanneal::LocalMinimizer::Minimize(bool dist_choice, int start_iters)
{
   int converged = 0;
   double new_energy;

   // first run
   energy = calc_new_energy();

   auto start = std::chrono::high_resolution_clock::now();

   std::pair<int,int> prev_pair(0,0);

   int iters = 1;

   doci2DM::BoundaryPoint *obj_bp = dynamic_cast<doci2DM::BoundaryPoint *> (method.get());

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

      int idx = 0;
      std::pair<int,int> tmp;

      if(dist_choice)
      {
         idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = 0;
      }

      tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

      // don't do the same pair twice in a row
      if(tmp==prev_pair)
         idx++;

      const auto& new_rot = list_rots[idx];
      prev_pair = std::make_pair(std::get<0>(new_rot), std::get<1>(new_rot));

      if(dist_choice)
         std::cout << iters << " (" << converged << ") Chosen: " << idx << std::endl;

      assert(ham->getOrbitalIrrep(std::get<0>(new_rot)) == ham->getOrbitalIrrep(std::get<1>(new_rot)));
      // do Jacobi rotation twice: once for the Hamiltonian data and once for the Unitary Matrix
      orbtrans->DoJacobiRotation(*ham, std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));
      orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));


      // start from zero every 25 iterations
      if(obj_bp && iters%25==0)
      {
         std::cout << "Restarting from zero" << std::endl;
         obj_bp->getX() = 0;
         obj_bp->getZ() = 0;
      }

      new_energy = calc_new_energy(*ham);

      std::stringstream h5_name;
      h5_name << getenv("SAVE_H5_PATH") << "/unitary-" << start_iters+iters << ".h5";
      orbtrans->get_unitary().saveU(h5_name.str());

      h5_name.str("");
      h5_name << getenv("SAVE_H5_PATH") << "/ham-" << start_iters+iters << ".h5";
      ham->save2(h5_name.str());

      h5_name.str("");
      h5_name << getenv("SAVE_H5_PATH") << "/rdm-" << start_iters+iters << ".h5";
      method->getRDM().WriteToFile(h5_name.str());

      if(obj_bp)
      {
         // if energy goes up instead of down, reset the start point
         if((std::get<3>(new_rot) - new_energy) < -1e-5)
         {
            std::cout << "Restarting from zero because too much up: " << std::get<3>(new_rot) - new_energy << std::endl;
            obj_bp->getX() = 0;
            obj_bp->getZ() = 0;
            obj_bp->Run();
            std::cout << "After restarting found: " << obj_bp->getEnergy() << " vs " << new_energy << "\t" << obj_bp->getEnergy()-new_energy << std::endl;
            new_energy = obj_bp->getEnergy();
         }

         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/X-" << start_iters+iters << ".h5";
         obj_bp->getX().WriteToFile(h5_name.str());

         h5_name.str("");
         h5_name << getenv("SAVE_H5_PATH") << "/Z-" << start_iters+iters << ".h5";
         obj_bp->getZ().WriteToFile(h5_name.str());
      }

      if(method->FullyConverged())
      {
         if(fabs(energy-new_energy)<conv_crit)
            converged++;
      }

      std::cout << iters << " (" << converged << ")\tRotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << std::get<3>(new_rot)+ham->getEconst() << "  E = " << new_energy+ham->getEconst() << "\t" << fabs(energy-new_energy) << std::endl;


      energy = new_energy;

      iters++;

      if(iters>1000)
      {
         std::cout << "Done 1000 steps, quiting..." << std::endl;
         break;
      }

      if(stopping_min)
         break;
   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Minimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   std::stringstream h5_name;
   h5_name << getenv("SAVE_H5_PATH") << "/optimale-uni.h5";
   get_Optimal_Unitary().saveU(h5_name.str());

   return iters;
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

/**
 * Choose a pair of orbitals to rotate over, according to the distribution of their relative
 * energy change.
 * @param orbs the list returned by scan_orbitals()
 * @return the index of the pair of orbitals in orbs
 */
int simanneal::LocalMinimizer::choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &orbs)
{
   std::uniform_real_distribution<double> dist(0, 1);

   const double choice = dist(mt);

   double norm = 0;

   for(auto &orb_pair: orbs)
      norm += (energy - std::get<3>(orb_pair));

   double cum = 0;
   for(int i=0;i<orbs.size();i++)
   {
      cum += (energy - std::get<3>(orbs[i]))/norm;
      if(choice < cum)
         return i;
   }

   assert(0 && "Should never ever be reached!");
   return -1;
}

int simanneal::LocalMinimizer::Minimize_noOpt(double stopcrit)
{
   int converged = 0;
   double old_energy;

   // first run
   energy = calc_new_energy();
   old_energy = energy;

   auto start = std::chrono::high_resolution_clock::now();

   std::pair<int,int> prev_pair(0,0);

   int iters = 1;

   while(fabs(old_energy-energy)<stopcrit && converged<conv_steps)
   {
      old_energy = energy;

      auto list_rots = scan_orbitals();

      std::sort(list_rots.begin(), list_rots.end(),
            [](const std::tuple<int,int,double,double> & a, const std::tuple<int,int,double,double> & b) -> bool
            {
            return std::get<3>(a) < std::get<3>(b);
            });

      for(auto& elem: list_rots)
         std::cout << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t" << std::get<3>(elem)+ham->getEconst() << "\t" << std::get<2>(elem) << std::endl;

      int idx = 0;
      std::pair<int,int> tmp;

      tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

      // don't do the same pair twice in a row
      if(tmp==prev_pair)
         idx++;

      const auto& new_rot = list_rots[idx];
      prev_pair = std::make_pair(std::get<0>(new_rot), std::get<1>(new_rot));

      assert(ham->getOrbitalIrrep(std::get<0>(new_rot)) == ham->getOrbitalIrrep(std::get<1>(new_rot)));
      // do Jacobi rotation twice: once for the Hamiltonian data and once for the Unitary Matrix
      orbtrans->DoJacobiRotation(*ham, std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));
      orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot));

      energy = std::get<3>(new_rot);

      std::cout << iters << " (" << converged << ")\tRotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << energy+ham->getEconst() << "\t" << old_energy-energy << std::endl;

      iters++;

      assert(energy<old_energy && "No minimization ?!?");

      if(fabs(old_energy-energy)<conv_crit)
         converged++;

      if(iters>10000)
      {
         std::cout << "Done 10000 steps, quiting..." << std::endl;
         break;
      }

      if(stopping_min)
         break;
   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Minimization with optimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   std::stringstream h5_name;
   h5_name << getenv("SAVE_H5_PATH") << "/optimale-uni-no-opt.h5";
   get_Optimal_Unitary().saveU(h5_name.str());

   return iters;
}

/**
 * Combine LocalMinimizer::Minimize and LocalMinimizer::Minimize_noOpt
 */
int simanneal::LocalMinimizer::Minimize_hybrid()
{
   const auto old_conv_crit = conv_crit;
   const auto old_conv_steps = conv_steps;
   const double switch_crit = 1e-1;
   int iters_noopt = 0;
   int iters = 0;

   while(iters_noopt < 20)
   {
      conv_crit = switch_crit;
      conv_steps = 10;

      iters += Minimize(false,iters);

      conv_crit = old_conv_crit;
      conv_steps = old_conv_steps;

      iters_noopt = Minimize_noOpt(switch_crit);
   }

   iters += Minimize(false,iters);

   return iters + iters_noopt;
}

/* vim: set ts=3 sw=3 expandtab :*/
