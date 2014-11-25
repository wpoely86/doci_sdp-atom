#include <iostream>
#include <getopt.h>

#include "include.h"

int main(int argc, char **argv)
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

   TPM ham(L,N);
   ham.HF_molecule(integralsfile);
   double norm_ham = std::sqrt(ham.ddot(ham));
   ham /= norm_ham;

   TPM rdm(L,N);
   rdm.init(lineq);

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

   lineq.check(rdm);

   int tot_iter = 0.0;

   double t = 1.0;
   double tolerance = 1.0e-5;
   double target = 1e-12;
   int iter = 0;

   TPM backup_rdm(rdm);

   //outer iteration: scaling of the potential barrier
   while(t > target)
   {
      cout << iter << "\t" << t << "\t" << rdm.getMatrices().trace() << "\t" << rdm.getVectors().trace() << "\t" << rdm.ddot(ham)*norm_ham + nuclrep << "\t" << rdm.S_2() << std::endl;

      double convergence = 1.0;
      iter++;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance)
      {
         tot_iter++;

         SUP P(L,N);

         P.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(L,N);

         grad.constr_grad(t,P,ham,lineq);

         //dit wordt de stap:
         TPM delta(L,N);

         //los het hessiaan stelsel op:
         cout << delta.solve(t,P,grad,lineq) << endl;

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);
      }

      cout << endl;
      t /= 2.0;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < target)
         tolerance = target;

      //extrapolatie:
      TPM extrapol(rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = rdm;

      double a = extrapol.line_search(t,rdm,ham);

      rdm.daxpy(a,extrapol);
   } 

   cout << endl;
   cout << "Energy: " << rdm.ddot(ham)*norm_ham + nuclrep << endl;
   cout << "Trace: " << rdm.trace() << endl;
   cout << "S^2: " << rdm.S_2() << endl;
   cout << "nuclrep: " << nuclrep << endl;

   cout << endl;
   cout << "total nr of iterations = " << tot_iter << endl;
   lineq.check(rdm);

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "rdm.h5";

   rdm.WriteToFile(h5_name.str().c_str());

   return 0;
}

/*  vim: set ts=3 sw=3 expandtab :*/
