#include <iostream>
#include <cmath>
#include <getopt.h>

// from CheMPS2
#include "Hamiltonian.h"
#include "OptIndex.h"


/**
 * This program generates the matrix elements for the
 * 1D Fermi-Hubbard model and saves them in a CheMPS2 
 * Hamiltonian class.
 */

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;

   std::string output;
   int L = 0;
   int N = 0;
   double U = 0;
   bool momspace = false;

   struct option long_options[] =
   {
      {"output",  required_argument, 0, 'o'},
      {"interaction",  required_argument, 0, 'U'},
      {"sites",  required_argument, 0, 'L'},
      {"particles",  required_argument, 0, 'N'},
      {"momspace",  no_argument, 0, 'm'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "ho:U:L:N:m", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -o, --output=output-filename    Set the output filename\n"
               "    -L, --sites=L                   Set the number of sites\n"
               "    -N, --particles=N               Set the number of particles\n"
               "    -U, --interaction=U             Set the on-site interaction\n"
               "    -m, --momspace                  Work in momentum space\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'o':
            output = optarg;
            break;
         case 'U':
            U = atof(optarg);
            break;
         case 'L':
            L = atoi(optarg);
            break;
         case 'N':
            N = atoi(optarg);
            break;
         case 'm':
            momspace = true;
            break;
      }

   if(! (L && N))
   {
      std::cerr << "You need to specifiy the system!" << endl;
      return 1;
   }

   cout << "Creating for L= " << L << " N= " << N << " U= " << U << endl;

   const std::vector<int> orb2irrep (L, 0);

   CheMPS2::Hamiltonian ham(L, 0, orb2irrep.data());
   // put everything to zero
   ham.reset();
   ham.setNe(N);
   ham.setEconst(0);

   if(momspace)
   {
      // one-particle integrals
      for(int i=0;i<L;i++)
         ham.setTmat(i, i, -2*std::cos(2*M_PI/(1.0*L)*i));

      // two-particle integrals
      for(int k1=0;k1<L;k1++)
         for(int k2=0;k2<L;k2++)
            for(int k3=0;k3<L;k3++)
               for(int k4=0;k4<L;k4++)
                  if((k1+k2)%L == (k3+k4)%L)
                     ham.setVmat(k1,k2,k3,k4, U*1.0/L);

   } else
   {
      // one-particle integrals
      for(int i=1;i<(L-1);i++)
      {
         ham.setTmat(i, i+1, -1);
         ham.setTmat(i, i-1, -1);
      }
      ham.setTmat(0, 1, -1);
      ham.setTmat(L-1, L-2, -1);

      // periodic boundary condition
      ham.setTmat(0, L-1, -1);
      ham.setTmat(L-1, 0, -1);

      // two-particle integrals
      for(int i=0;i<L;i++)
         ham.setVmat(i, i, i, i, U);
   }

   if(output.empty())
   {
      output = "hub-integrals-";
      if(momspace)
         output += "mom-";
      
      output += std::to_string(L) + "-" + std::to_string(N) + "-" + std::to_string(U) + ".h5";
   }

   cout << "Writing Hamiltonian to " << output << endl;

   ham.save2(output);

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
