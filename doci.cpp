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

#include <iostream>
#include <chrono>
#include <getopt.h>
#include <mpi.h>
#include <signal.h>

#include "include.h"
#include "BoundaryPoint.h"
#include "PotentialReducation.h"
// from CheMPS2
#include "Hamiltonian.h"

#include "SimulatedAnnealing.h"

sig_atomic_t stopping = 0;
sig_atomic_t stopping_min = 0;

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;
   using namespace doci2DM;
   using simanneal::SimulatedAnnealing;

   cout.precision(10);

   std::string integralsfile = "mo-integrals.h5";
   bool bp = false;
   bool pr = false;
   std::string startinput;

   struct option long_options[] =
   {
      {"integrals",  required_argument, 0, 'i'},
      {"start",  required_argument, 0, 's'},
      {"boundary-point",  no_argument, 0, 'b'},
      {"potential-reduction",  no_argument, 0, 'p'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "hi:bps:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -i, --integrals=integrals-file  Set the input integrals file\n"
               "    -b, --boundary-point            Use the boundary point method as solver (default)\n"
               "    -p, --potential-reduction       Use the potential reduction method as solver\n"
               "    -s, --start                     Use this a start point for the Simulated Annealing\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'i':
            integralsfile = optarg;
            break;
         case 's':
            startinput = optarg;
            break;
         case 'b':
            bp = true;
            break;
         case 'p':
            pr = true;
            break;
      }

   if(!bp && !pr)
   {
      std::cout << "Tell me what to do..." << std::endl;
      return 0;
   }

   cout << "Reading: " << integralsfile << endl;

   // make sure we have a save path, even if it's not specify already
   // This will not overwrite an already set SAVE_H5_PATH
   setenv("SAVE_H5_PATH", "./", 0);

   cout << "Using save path: " << getenv("SAVE_H5_PATH") << endl;

   MPI_Init(&argc,&argv);

   SimulatedAnnealing opt(CheMPS2::Hamiltonian::CreateFromH5(integralsfile));

   if(bp)
      opt.UseBoundaryPoint();
   else if(pr)
      opt.UsePotentialReduction();

   auto& ham = opt.getHam();

   const auto L = ham.getL(); //dim sp hilbert space
   const auto N = ham.getNe(); //nr of particles

   cout << "Starting with L=" << L << " N=" << N << endl;
   auto start = std::chrono::high_resolution_clock::now();

   opt.Set_start_temp(0.1);
   opt.Set_delta_temp(0.99);
   opt.Set_max_angle(1.3);
   opt.Set_delta_angle(0.999);

//   opt.get_Optimal_Unitary().loadU("optimale-uni.h5");
//   opt.calc_new_energy();

   if(bp)
   {
      try {
         opt.getMethod_BP().set_use_prev_result(true);
      } catch (const std::bad_cast err)
      {
         std::cerr << "Shit just hit the fan!\t" << err.what() << std::endl;
      }
   }

   if(! startinput.empty())
      opt.getMethod().getRDM().ReadFromFile(startinput.c_str());
   else
      opt.calc_energy();

   opt.optimize_mpi();

   auto end = std::chrono::high_resolution_clock::now();

   cout << "Total Runtime: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << endl;

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if(rank == 0)
   {
      std::stringstream h5_name1;
      h5_name1 << getenv("SAVE_H5_PATH") << "/rdm.h5";
      opt.getMethod().getRDM().WriteToFile(h5_name1.str().c_str());

      std::stringstream h5_name2;
      h5_name2 << getenv("SAVE_H5_PATH") << "/optimale-uni.h5";
      opt.get_Optimal_Unitary().saveU(h5_name2.str());
   }

   MPI_Finalize();

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
