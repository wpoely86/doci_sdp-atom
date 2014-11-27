#include <iostream>
#include <getopt.h>

#include "include.h"
#include "BoundaryPoint.h"
#include "PotentialReducation.h"
// from CheMPS2
#include "Hamiltonian.h"

#include "SimulatedAnnealing.h"

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

   struct option long_options[] =
   {
      {"integrals",  required_argument, 0, 'i'},
      {"boundary-point",  no_argument, 0, 'b'},
      {"potential-reduction",  no_argument, 0, 'p'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "hi:bp", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -i, --integrals=integrals-file  Set the input integrals file\n"
               "    -b, --boundary-point            Use the boundary point method as solver (default)\n"
               "    -p, --potential-reduction       Use the potential reduction method as solver\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'i':
            integralsfile = optarg;
            break;
         case 'b':
            bp = true;
            break;
         case 'p':
            pr = true;
            break;
      }

   cout << "Reading: " << integralsfile << endl;

   SimulatedAnnealing opt(CheMPS2::Hamiltonian::CreateFromH5(integralsfile));

   if(bp)
      opt.UseBoundaryPoint();
   else if(pr)
      opt.UsePotentialReduction();

   auto& ham = opt.getHam();

   const auto L = ham.getL(); //dim sp hilbert space
   const auto N = ham.getNe(); //nr of particles

   cout << "Starting with L=" << L << " N=" << N << endl;

   opt.Set_start_temp(0.1);
   opt.Set_delta_temp(0.99);
   opt.Set_max_angle(1.3);
   opt.Set_delta_angle(0.999);

//   opt.get_Optimal_Unitary().loadU("optimale-uni.h5");
//   opt.calc_new_energy();

   opt.calc_energy();

   opt.optimize();

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "rdm.h5";

   opt.getMethod().getRDM().WriteToFile(h5_name.str().c_str());

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
