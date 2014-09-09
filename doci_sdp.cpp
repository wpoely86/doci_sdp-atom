#include <iostream>
#include <getopt.h>

#include "include.h"
#include "integrals.h"

int main(int argc, char **argv)
{
   using std::cout;
   using std::endl;

   cout.precision(10);

   std::string inputfile = "start.stp";

   std::string integralsfile = "mo-integrals.h5";

   struct option long_options[] =
   {
      {"file",  required_argument, 0, 'f'},
      {"integrals",  required_argument, 0, 'i'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "hf:i:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -f, --file=input-file           Set the input file\n"
               "    -i, --integrals=integrals-file  Set the input integrals file\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'f':
            inputfile = optarg;
            break;
         case 'i':
            integralsfile = optarg;
            break;
      }

   cout << "Reading: " << inputfile << endl;

   integrals::init(inputfile.c_str());

   //here the cartesian integrals are calculated
   integrals::CartInt::calc_integrals();

   integrals::CartInt::orthogonalize();

   const int L = integrals::CI_SPM::gdim();//dim sp hilbert space
   const int N = integrals::input::NumberOfElectrons();//nr of particles

   cout << "Starting with L=" << L << " N=" << N << endl;





   integrals::clean();

   return 0;
}

/*  vim: set ts=3 sw=3 expandtab :*/
