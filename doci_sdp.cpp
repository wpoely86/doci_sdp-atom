#include <iostream>
#include <getopt.h>

#include "include.h"
#include "PotentialReducation.h"
#include "Hamiltonian.h"

int main(int argc, char **argv)
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

   auto ham = CheMPS2::Hamiltonian::CreateFromH5(integralsfile);
   PotentialReduction method(ham);

   const auto L = ham.getL(); //nr of particles
   const auto N = ham.getNe(); //nr of particles

   cout << "Starting with L=" << L << " N=" << N << endl;

//   Matrix rdmfull(rdm.gn());
//   rdmfull = 0;
//
//   for(int i=0;i<L;i++)
//      for(int j=0;j<L;j++)
//         rdmfull(i,j) = rdm.getMatrix(0)(i,j);
//
//   int tmp_dim = rdm.getVector(0).gn();
//   for(int a=0;a<4;a++)
//      for(int i=0;i<tmp_dim;i++)
//         rdmfull(L+a*tmp_dim+i,L+a*tmp_dim+i) = rdm.getVector(0)[i];
//   rdmfull.SaveRawToFile("rdm-OK-full-P.h5");

   method.Run();

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "rdm.h5";

   method.getRDM().WriteToFile(h5_name.str().c_str());

   return 0;
}

/*  vim: set ts=3 sw=3 expandtab :*/
