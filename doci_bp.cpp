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
#include <fstream>
#include <cstring>
#include <getopt.h>
#include <signal.h>

#include "include.h"
#include "BoundaryPoint.h"
#include "LocalMinimizer.h"

// from CheMPS2
#include "Hamiltonian.h"
#include "OptIndex.h"

// if set, the signal has been given to stop the calculation and write current step to file
sig_atomic_t stopping = 0;
sig_atomic_t stopping_min = 0;

void stopcalcsignal(int sig);
void stopminsignal(int sig);

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;
   using namespace doci2DM;
   using simanneal::LocalMinimizer;

   cout.precision(10);

   std::string integralsfile = "mo-integrals.h5";
   std::string unitary;
   std::string rdmfile;
   bool random = false;
   bool localmini = false;
   bool scan = false;
   bool localmininoopt = false;

   struct option long_options[] =
   {
      {"integrals",  required_argument, 0, 'i'},
      {"unitary",  required_argument, 0, 'u'},
      {"rdm",  required_argument, 0, 'd'},
      {"random",  no_argument, 0, 'r'},
      {"scan",  no_argument, 0, 's'},
      {"local-minimizer",  no_argument, 0, 'l'},
      {"local-minimizer-no-opt",  no_argument, 0, 'n'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "d:rlhi:u:sn", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -i, --integrals=integrals-file  Set the input integrals file\n"
               "    -d, --rdm=rdm-file              Use this rdm as starting point\n"
               "    -u, --unitary=unitary-file      Use the unitary matrix in this file\n"
               "    -r, --random                    Perform a random unitary transformation on the Hamiltonian\n"
               "    -l, --local-minimizer           Use the local minimizer\n"
               "    -n, --local-minimizer-no-opt    Use the local minimizer without optimalization\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'i':
            integralsfile = optarg;
            break;
         case 'u':
            unitary = optarg;
            break;
         case 'd':
            rdmfile = optarg;
            break;
         case 'r':
            random = true;
            break;
         case 'l':
            localmini = true;
            break;
         case 's':
            scan = true;
            break;
         case 'n':
            localmininoopt = true;
            break;
      }

   cout << "Reading: " << integralsfile << endl;

   // make sure we have a save path, even if it's not specify already
   // This will not overwrite an already set SAVE_H5_PATH
   setenv("SAVE_H5_PATH", "./", 0);

   cout << "Using save path: " << getenv("SAVE_H5_PATH") << endl;

   auto ham = CheMPS2::Hamiltonian::CreateFromH5(integralsfile);

   const auto L = ham.getL(); //dim sp hilbert space
   const auto N = ham.getNe(); //nr of particles

   cout << "Starting with L=" << L << " N=" << N << endl;

   if(!unitary.empty() && !localmini)
   {
      cout << "Reading transform: " << unitary << endl;

      simanneal::OrbitalTransform orbtrans(ham);

      orbtrans.get_unitary().loadU(unitary);
      orbtrans.fillHamCI(ham);
   }

   const simanneal::OptIndex opt(ham);

   if(random)
   {
      simanneal::UnitaryMatrix X(opt);
      X.fill_random();
      X.make_skew_symmetric();
      X.print_unitary();
      simanneal::OrbitalTransform orbtrans(ham);
      orbtrans.update_unitary(X, false);
      std::string filename = getenv("SAVE_H5_PATH");
      filename += "/random-start-unitary.h5";
      orbtrans.get_unitary().saveU(filename);
      orbtrans.fillHamCI(ham);
   }


   BoundaryPoint method(ham);
   method.set_tol_PD(1e-7);
//   method.getLineq() = Lineq(L,N,true);

   char *X_env = getenv("v2DM_DOCI_SUP_X");
   if(X_env && strlen(X_env) > 0)
   {
      cout << "Reading X from " << X_env << endl;
      method.getX().ReadFromFile(X_env);
   }

   char *Z_env = getenv("v2DM_DOCI_SUP_Z");
   if(Z_env && strlen(Z_env) > 0)
   {
      cout << "Reading Z from " << Z_env << endl;
      method.getZ().ReadFromFile(Z_env);
   }

   // set up everything to handle SIGALRM
   struct sigaction act;
   act.sa_flags = 0;
   act.sa_handler = &stopcalcsignal;

   sigset_t blockset;
   sigemptyset(&blockset); // we don't block anything in the handler
   act.sa_mask = blockset;

   sigaction(SIGALRM, &act, 0);

   act.sa_handler = &stopminsignal;

   sigaction(SIGUSR1, &act, 0);

   if(!rdmfile.empty())
   {
      cout << "Reading rdm: " << rdmfile << endl;
      method.getRDM().ReadFromFile(rdmfile);
   }

   if(localmini)
   {
      LocalMinimizer minimize(ham);

      if(!unitary.empty())
      {
         cout << "Starting local minimizer from: " << unitary << endl;
         minimize.getOrbitalTf().get_unitary().loadU(unitary);
      }

      minimize.UseBoundaryPoint();
      minimize.getMethod_BP().getX() = method.getX();
      minimize.getMethod_BP().getZ() = method.getZ();
      minimize.getMethod_BP().set_use_prev_result(true);
      minimize.getMethod_BP().set_tol_PD(1e-7);
//      minimize.getMethod_BP().getLineq() = Lineq(L,N,true);
      minimize.set_conv_steps(10);
//      minimize.getMethod_BP().set_max_iter(5);

      minimize.set_conv_crit(1e-6);

      if(localmininoopt)
         minimize.Minimize_noOpt(1e-2);
      else
         minimize.Minimize();

      cout << "Bottom is " << minimize.get_energy() << endl;

      method = minimize.getMethod_BP();
      ham = minimize.getHam();
   }

   method.set_use_prev_result(false);
   method.Reset_avg_iters();
   method.Run();

   cout << "The optimal energy is " << method.evalEnergy() << std::endl;

   if(scan)
      Tools::scan_all_bp(method.getRDM(), ham);



/* //   for(int k_in=0;k_in<L;k_in++)
 * //      for(int l_in=k_in+1;l_in<L;l_in++)
 *    int k_in = 2;
 *    int l_in = 3;
 *          if(ham.getOrbitalIrrep(k_in) == ham.getOrbitalIrrep(l_in))
 *          {
 *             std::fstream fs;
 *             std::string filename = "orbs-scan-" + std::to_string(k_in) + "-" + std::to_string(l_in) + ".txt";
 *             fs.open(filename, std::fstream::out | std::fstream::trunc);
 * 
 *             fs.precision(10);
 * 
 *             fs << "# theta\trot\trot+v2dm" << endl;
 * 
 *             auto found = rdm.find_min_angle(k_in,l_in,0.3,getT3,getV3);
 * 
 *             cout << "Min:\t" << k_in << "\t" << l_in << "\t" << found.first << "\t" << found.second << endl;
 * 
 *             cout << "######################" << endl;
 * 
 *             int Na = 200;
 * #pragma omp parallel for
 *             for(int a=0;a<=Na;a++)
 *             {
 *                double theta = 2.0*M_PI/(1.0*Na) * a;
 * 
 *                BoundaryPoint mymethod(ham);
 * 
 *                mymethod.getRDM() = orig_rdm;
 *                mymethod.getHam() = orig_ham;
 *                mymethod.getHam().rotate(k_in, l_in, theta, getT3, getV3);
 * 
 *                simanneal::OrbitalTransform orbtrans2(ham2);
 *                orbtrans2.DoJacobiRotation(ham, k_in, l_in, theta);
 *                orbtrans2.get_unitary().jacobi_rotation(ham2.getOrbitalIrrep(k_in), k_in, l_in, theta);
 *                //orbtrans2.fillHamCI(ham);
 * 
 *                mymethod.BuildHam(new_ham);
 * 
 *                double new_en = orig_rdm.ddot(mymethod.getHam());
 * 
 *                mymethod.Run();
 * 
 * #pragma omp critical
 *                fs << theta << "\t" << new_en << "\t" << mymethod.getEnergy() << endl;
 *             }
 *             fs.close();
 *          }
 */


   std::string h5_name = getenv("SAVE_H5_PATH");
   h5_name += "/optimal-rdm.h5";

   method.getRDM().WriteToFile(h5_name);

   h5_name = getenv("SAVE_H5_PATH");
   h5_name += "/optimal-ham.h5";

   ham.save2(h5_name);

   return 0;
}

void stopcalcsignal(int sig)
{
   stopping=1;
}

void stopminsignal(int sig)
{
   stopping_min=1;
}

/* vim: set ts=3 sw=3 expandtab :*/
