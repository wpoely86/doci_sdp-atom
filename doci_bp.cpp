#include <iostream>
#include <fstream>
#include <getopt.h>

#include "include.h"
#include "BoundaryPoint.h"
// from CheMPS2
#include "Hamiltonian.h"
#include "SimulatedAnnealing.h"
#include "OptIndex.h"
#include "LocalMinimizer.h"

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;
   using namespace doci2DM;
   using simanneal::SimulatedAnnealing;
   using simanneal::LocalMinimizer;

   cout.precision(10);

   std::string integralsfile = "mo-integrals.h5";
   std::string unitary;
   std::string rdmfile;
   bool random = false;
   bool localmini = false;

   struct option long_options[] =
   {
      {"integrals",  required_argument, 0, 'i'},
      {"unitary",  required_argument, 0, 'u'},
      {"rdm",  required_argument, 0, 'd'},
      {"random",  no_argument, 0, 'r'},
      {"local-minimizer",  no_argument, 0, 'l'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "d:rlhi:u:", long_options, &i)) != -1)
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

   if(! unitary.empty())
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
   auto &rdm = method.getRDM();

   if(!rdmfile.empty())
   {
      cout << "Reading rdm: " << rdmfile << endl;
      method.getRDM().ReadFromFile(rdmfile);
   } else
   method.Run();

   if(localmini)
   {
      LocalMinimizer minimize(ham);

      minimize.getMethod_BP() = method;

      minimize.getMethod_BP().set_use_prev_result(true);
      minimize.getMethod_BP().set_tol_PD(1e-7);
      minimize.set_conv_crit(1e-5);

      minimize.Minimize();

      cout << "Bottom is " << minimize.get_energy() << endl;

      method.getRDM() = minimize.getMethod().getRDM();
      method.getHam() = minimize.getMethod().getHam();
      ham = minimize.getHam();
   }


//   method.set_use_prev_result(true);
   method.Run();

   cout << "The optimal energy is " << method.evalEnergy() << std::endl;

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

/* vim: set ts=3 sw=3 expandtab :*/

