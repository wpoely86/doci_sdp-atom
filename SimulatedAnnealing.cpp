#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <unistd.h>
#include <hdf5.h>
#include <mpi.h>

#include "SimulatedAnnealing.h"
#include "OptIndex.h"
#include "BoundaryPoint.h"
#include "PotentialReducation.h"

/**
 * You still need to set the max_angle, delta_angle, start_temp and
 * delta_temp after creating the object.
 * @param mol the molecular data to use
 */
simanneal::SimulatedAnnealing::SimulatedAnnealing(const CheMPS2::Hamiltonian &mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   OptIndex index(*ham);

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   mt = std::mt19937_64(rd());

   steps = 0;
   energy = 0;
   max_steps = 20000;
   unaccepted = 0;
}

simanneal::SimulatedAnnealing::SimulatedAnnealing(CheMPS2::Hamiltonian &&mol)
{
   ham.reset(new CheMPS2::Hamiltonian(mol));

   OptIndex index(*ham);

   opt_unitary.reset(new UnitaryMatrix(index));

   orbtrans.reset(new OrbitalTransform(*ham));

   method.reset(new doci2DM::BoundaryPoint(*ham));

   mt = std::mt19937_64(rd());

   steps = 0;
   energy = 0;
   max_steps = 20000;
   unaccepted = 0;
}

simanneal::SimulatedAnnealing::~SimulatedAnnealing() = default;

/**
 * @return the real energy (calculated + nuclear repulsion)
 */
double simanneal::SimulatedAnnealing::get_energy() const
{
   return energy + ham->getEconst();
}

void simanneal::SimulatedAnnealing::Set_max_angle(double max_angle)
{
   this->max_angle = max_angle;
}

void simanneal::SimulatedAnnealing::Set_delta_angle(double delta_angle)
{
   this->delta_angle = delta_angle;
}

void simanneal::SimulatedAnnealing::Set_start_temp(double start_temp)
{
   this->start_temp = start_temp;
}

void simanneal::SimulatedAnnealing::Set_delta_temp(double delta_temp)
{
   this->delta_temp = delta_temp;
}

/**
 * Decide wether or not to accept the new energy
 * @param e_new the new energy
 * @return accept or not
 */
bool simanneal::SimulatedAnnealing::accept_function(double e_new)
{
   std::uniform_real_distribution<double> dist_accept(0, 1);

   if(e_new < energy)
      return true;
   else
   {
      double chance = std::exp((energy - e_new) / cur_temp);

      if ( dist_accept(mt) * (1+chance) > chance)
         //biggest chance for failing
         return false;
      else 
         return true;
   }
}

/**
 * Calculate the energy with the current
 * molecular data
 */
void simanneal::SimulatedAnnealing::calc_energy()
{
   method->BuildHam(*ham);
   method->Run();

   energy = method->getEnergy();
}

/**
 * Calculate the energy with the current
 * molecular data
 */
double simanneal::SimulatedAnnealing::calc_new_energy()
{
   orbtrans->fillHamCI(*ham);

   method->BuildHam(*ham);
   method->Run();

   return method->getEnergy();
}

/**
 * Do the simulated annealing
 */
void simanneal::SimulatedAnnealing::optimize()
{
   std::uniform_int_distribution<int> dist(0, ham->getL()-1);
   std::uniform_real_distribution<double> dist_angles(0, 1);

   unaccepted = 0;
   double lowest_energy = energy;
   cur_temp = start_temp;

   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   std::stringstream out_name;
   out_name << getenv("SAVE_H5_PATH") << "/output-" << rank << ".txt";

   std::ostream* fp = &std::cout;
   std::ofstream fout;
   if(size > 1)
   {
      fout.open(out_name.str(), std::ios::out | std::ios::app);
      fp = &fout;
   }
   std::ostream &out = *fp;
   out.precision(10);

   helpers::matrix sample_pairs(ham->getL(), ham->getL());
   sample_pairs = 0;

   for(int i=0;i<ham->getL();i++)
      for(int j=i+1;j<ham->getL();j++)
         if(ham->getOrbitalIrrep(i) != ham->getOrbitalIrrep(j))
            sample_pairs(i,j) = sample_pairs(j,i) = -1;

   for(unsigned int i=0;i<max_steps;i++)
   {
      auto orb1 = dist(mt);
      auto orb2 = dist(mt);

      if(orb1 != orb2 && ham->getOrbitalIrrep(orb1) == ham->getOrbitalIrrep(orb2))
      {
         sample_pairs(orb1, orb2) += 1;
         sample_pairs(orb2, orb1) += 1;

         // between -1 and 1 but higher probablity to be close to zero (seems to work better)
         double cur_angle = max_angle * (dist_angles(mt) - dist_angles(mt));

         out << "P=" << rank << "\t" << i << "\tT=" << cur_temp << "\tOrb1=" << orb1 << "\tOrb2=" << orb2 << "  Over " << cur_angle << std::endl;

         orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(orb1), orb1, orb2, cur_angle);

         auto new_energy = calc_new_energy();

         if(new_energy < lowest_energy)
            lowest_energy = new_energy;

         out << "P=" << rank << "\tT=" << cur_temp << "\tNew energy = " << new_energy + ham->getEconst() << "\t Old energy = " << get_energy();

         if(accept_function(new_energy))
         {
            energy = new_energy;
            out << "\t=> Accepted" << std::endl;
         }
         else
         {
            unaccepted++;
            out << "\t=> Unaccepted, " << unaccepted << std::endl;
            orbtrans->get_unitary().jacobi_rotation(ham->getOrbitalIrrep(orb1), orb1, orb2, -1*cur_angle);
         }

         cur_temp *= delta_temp;
         max_angle *= delta_angle;

         if(unaccepted > 1500)
         {
            out << "Too many unaccepted, stopping" << std::endl;
            break;
         }

      }

   }

   out << "Bottom was " << lowest_energy + ham->getEconst() << std::endl;
   out << "Final energy = " << get_energy() << std::endl;

   for(int i=0;i<ham->getL();i++)
      for(int j=i+1;j<ham->getL();j++)
         if(sample_pairs(i,j) >= 0)
            out << ham->getOrbitalIrrep(i) << "\t" << i << "\t" << j << "\t" << sample_pairs(i,j) << std::endl;


   if(size > 1)
      fout.close();
}

simanneal::UnitaryMatrix& simanneal::SimulatedAnnealing::get_Optimal_Unitary()
{
   return orbtrans->get_unitary();
}

CheMPS2::Hamiltonian& simanneal::SimulatedAnnealing::getHam() const
{
   return *ham;
}

simanneal::OrbitalTransform& simanneal::SimulatedAnnealing::getOrbitalTf() const
{
   return *orbtrans;
}

doci2DM::Method& simanneal::SimulatedAnnealing::getMethod() const
{
   return *method;
}

void simanneal::SimulatedAnnealing::UseBoundaryPoint()
{
   method.reset(new doci2DM::BoundaryPoint(*ham));
}

void simanneal::SimulatedAnnealing::UsePotentialReduction()
{
   method.reset(new doci2DM::PotentialReduction(*ham));
}

void simanneal::SimulatedAnnealing::optimize_mpi()
{
   max_steps = 200;
   bool keepsearch = true;

   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   std::stringstream out_name;
   out_name << getenv("SAVE_H5_PATH") << "/output-" << rank << ".txt";
   if(size > 1)
      method->set_outfile(out_name.str());
   std::ofstream myfile;
   myfile.open(out_name.str(), std::ios::out | std::ios::trunc | std::ios::binary);
   myfile.close();

   std::vector<double> energies(size);

   while(keepsearch)
   {
      auto start = std::chrono::high_resolution_clock::now();

      optimize();

      auto end = std::chrono::high_resolution_clock::now();

      MPI_Gather(&energy, 1, MPI_DOUBLE, energies.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      int min_rank = 0;

      if(rank == 0)
      {
         for(int i=0;i<size;i++)
         {
            std::cout << i << "\t" << energies[i] + ham->getEconst() << std::endl;
            if(energies[i] < energies[min_rank])
               min_rank = i;
         }

         energy = energies[min_rank];
         std::cout << "Rank " << min_rank << " has the lowest energy => " << get_energy() << std::endl;
      }

      MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&min_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

      std::stringstream h5_name;
      h5_name << getenv("SAVE_H5_PATH") << "/unitary-" << max_steps << "-" << rank << ".h5";
      orbtrans->get_unitary().saveU(h5_name.str());

      orbtrans->get_unitary().sendreceive(min_rank);

      max_steps *= 10;

      if(rank == 0)
         std::cout << "New max steps is " << max_steps << std::endl;

      start_temp = cur_temp / delta_temp;

      if(max_steps >= 200000)
         keepsearch = false;

      // sleep rank seconds to avoid overlap in printing (nasty, I know)
      sleep(rank+1);
      std::cout << "P=" << rank << "\tMPI Step runtime: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;
   }
}

/* vim: set ts=3 sw=3 expandtab :*/
