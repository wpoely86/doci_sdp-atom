#include <iostream>
#include <getopt.h>

#include "include.h"

int main(int argc,char **argv)
{
   using std::cout;
   using std::endl;

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

   const int L = Tools::getspDimension(integralsfile);//dim sp hilbert space
   const int N = Tools::getNumberOfParticles(integralsfile);//nr of particles
   const double nuclrep = Tools::getNuclearRepulEnergy(integralsfile);

   cout << "Starting with L=" << L << " N=" << N << endl;

   Lineq lineq(L,N);

   //hamiltoniaan
   TPM ham(L,N);
   ham.HF_molecule(integralsfile);

   TPM ham_copy(ham);

   //only traceless hamiltonian needed in program.
   ham.Proj_E(lineq);

   //primal
   SUP X(L,N);

   //dual
   SUP Z(L,N);

   //Lagrange multiplier
   SUP V(L,N);

   //just dubya
   SUP W(L,N);

   SUP u_0(L,N);

   //little help
   TPM hulp(L,N);

   u_0.init_S(lineq);

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tol_PD = 1.0e-6;
   double tol_en = 1.0e-3;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   // mazziotti uses 1.6 for this
   double mazzy = 1.0;

   int iter_dual,iter_primal(0);
   int max_iter = 2;

   int tot_iter = 0;

   while(P_conv > tol_PD || D_conv > tol_PD || fabs(convergence) > tol_en)
   {
      ++iter_primal;

      D_conv = 1.0;

      iter_dual = 0;

      while(D_conv > tol_PD  && iter_dual <= max_iter)
      {
         ++tot_iter;

         ++iter_dual;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(1.0/sigma,X);

         TPM b(L,N);

         b.collaps(B, lineq);

         b.daxpy(-1.0/sigma,ham);

         hulp.S(b, true);
//         hulp.InverseS(b, lineq);
//
         hulp.Proj_E(lineq);

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-mazzy/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v(L,N);

         v.collaps(V, lineq);

         v -= ham;

         D_conv = sqrt(v.ddot(v));
     }

      //update primal:
      X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z.getI().ddot(ham) + X.ddot(u_0);

      cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z.getI().ddot(ham_copy) + nuclrep << "\t" << Z.getI().S_2() << endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;
   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.getI()) + nuclrep << endl;
   cout << "Trace: " << Z.getI().trace() << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "S^2: " << Z.getI().S_2() << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   cout << endl;
   cout << "total nr of iterations = " << tot_iter << endl;

   lineq.check(Z.getI());

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "rdm.h5";

   Z.getI().WriteToFile(h5_name.str().c_str());

   return 0;
}

/* vim: set ts=3 sw=3 expandtab :*/
