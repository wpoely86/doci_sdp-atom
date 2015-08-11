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
   double g = 0;

   struct option long_options[] =
   {
      {"output",  required_argument, 0, 'o'},
      {"pairing",  required_argument, 0, 'g'},
      {"sites",  required_argument, 0, 'L'},
      {"particles",  required_argument, 0, 'N'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;

   while( (j = getopt_long (argc, argv, "ho:g:L:N:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -o, --output=output-filename    Set the output filename\n"
               "    -L, --sites=L                   Set the number of sites\n"
               "    -N, --particles=N               Set the number of particles\n"
               "    -g, --pairing=g                 Set the pairing strength\n"
               "    -h, --help                      Display this help\n"
               "\n";
            return 0;
            break;
         case 'o':
            output = optarg;
            break;
         case 'g':
            g = atof(optarg);
            break;
         case 'L':
            L = atoi(optarg);
            break;
         case 'N':
            N = atoi(optarg);
            break;
      }

   if(! (L && N))
   {
      std::cerr << "You need to specifiy the system!" << endl;
      return 1;
   }

   cout << "Creating for L= " << L << " N= " << N << " g= " << g << endl;

   const std::vector<int> orb2irrep (L, 0);

   CheMPS2::Hamiltonian ham(L, 0, orb2irrep.data());
   // put everything to zero
   ham.reset();
   ham.setNe(N);
   ham.setEconst(0);

   // one-particle integrals
   for(int i=0;i<L;i++)
      ham.setTmat(i, i, i+1);

   // two-particle integrals
   for(int i=0;i<L;i++)
      for(int j=0;j<L;j++)
         ham.setVmat(i, i, j, j, g);

   if(output.empty())
      output = "pairing-integrals-" + std::to_string(L) + "-" + std::to_string(N) + "-" + std::to_string(g) + ".h5";

   cout << "Writing Hamiltonian to " << output << endl;

   ham.save2(output);

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
