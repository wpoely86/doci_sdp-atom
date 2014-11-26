#include <iostream>
#include <getopt.h>

#include "include.h"
#include "BoundaryPoint.h"
// from CheMPS2
#include "Hamiltonian.h"

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;
   using namespace doci2DM;

   cout.precision(10);

   std::string integralsfile = "mo-integrals.h5";

   struct option long_options[] =
   {
      {"integrals",  required_argument, 0, 'i'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "hi:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -i, --integrals=integrals-file  Set the input integrals file\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'i':
            integralsfile = optarg;
            break;
      }

   cout << "Reading: " << integralsfile << endl;

//   const int L = Tools::getspDimension(integralsfile);//dim sp hilbert space
   const int N = Tools::getNumberOfParticles("CheMPS2_Ham_parent.h5");//nr of particles

   CheMPS2::Hamiltonian ham(true, "CheMPS2_Ham_parent.h5", "CheMPS2_Ham_Tmat.h5", "CheMPS2_Ham_Vmat.h5");

   cout << "Starting with L=" << ham.getL() << " N=" << N << endl;

   BoundaryPoint method(ham);

   method.Run();

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "rdm.h5";

   method.getZ().getI().WriteToFile(h5_name.str().c_str());

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
