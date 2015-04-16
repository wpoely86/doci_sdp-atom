#include <fstream>
#include <hdf5.h>
#include "include.h"

#include "PotentialReducation.h"
#include "BoundaryPoint.h"
#include "Hamiltonian.h"

using namespace doci2DM;

int Tools::getspDimension(std::string filename)
{
   hid_t       file_id, group_id, attribute_id;
   herr_t      status;
   int spdim = 0; // to avoid uninitalized errors

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/integrals", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Aopen(group_id, "sp_dim", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_INT, &spdim);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   return spdim;
}

int Tools::getNumberOfParticles(std::string filename)
{
   hid_t       file_id, group_id, attribute_id;
   herr_t      status;
   int nelectrons = 0; // to avoid uninitalized errors

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/Data", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Dopen(group_id, "nelectrons", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Dread(attribute_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nelectrons);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   return nelectrons;
}

double Tools::getNuclearRepulEnergy(std::string filename)
{
   hid_t       file_id, group_id, attribute_id;
   herr_t      status;
   double nuclrep = 0; // to avoid uninitalized errors

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/integrals", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Aopen(group_id, "nuclear_repulsion_energy", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &nuclrep);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   return nuclrep;
}


void Tools::scan_all(const TPM &rdm, const CheMPS2::Hamiltonian &ham)
{
   const int L = rdm.gL();

   std::function<double(int,int)> getT = [&ham] (int a, int b) -> double { return ham.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&ham]  (int a, int b, int c, int d) -> double { return ham.getVmat(a,b,c,d); };

   PotentialReduction mymethod(ham);
   
   auto orig_ham = mymethod.getHam();

   for(int k_in=0;k_in<L;k_in++)
      for(int l_in=k_in+1;l_in<L;l_in++)
         if(ham.getOrbitalIrrep(k_in) == ham.getOrbitalIrrep(l_in))
         {
            std::fstream fs;
            std::string filename = getenv("SAVE_H5_PATH");
            filename += "/orbs-scan-" + std::to_string(k_in) + "-" + std::to_string(l_in) + ".txt";
            fs.open(filename, std::fstream::out | std::fstream::trunc);

            fs.precision(10);

            fs << "# theta\trot\trot+v2dm" << std::endl;

            auto found = rdm.find_min_angle(k_in,l_in,0.3,getT,getV);

            std::cout << "Min:\t" << k_in << "\t" << l_in << "\t" << found.first << "\t" << found.second << std::endl;

            std::cout << "######################" << std::endl;


            int Na = 100;
//#pragma omp parallel for
            for(int a=0;a<=Na;a++)
            {
               double theta = 1.0*M_PI/(1.0*Na) * a - M_PI/2.0;

               mymethod.getHam() = orig_ham;
               mymethod.getRDM() = rdm;
               mymethod.getHam().rotate(k_in, l_in, theta, getT, getV);

               double new_en = mymethod.evalEnergy();

               mymethod.Run();

//#pragma omp critical
               fs << theta << "\t" << new_en << "\t" << mymethod.evalEnergy() << std::endl;
            }

            fs.close();
         }
}

void Tools::scan_all_bp(const TPM &rdm, const CheMPS2::Hamiltonian &ham)
{
   const int L = rdm.gL();

   std::function<double(int,int)> getT = [&ham] (int a, int b) -> double { return ham.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&ham]  (int a, int b, int c, int d) -> double { return ham.getVmat(a,b,c,d); };

   BoundaryPoint method(ham);

   const auto orig_ham = method.getHam();

   for(int k_in=0;k_in<L;k_in++)
      for(int l_in=k_in+1;l_in<L;l_in++)
         if(ham.getOrbitalIrrep(k_in) == ham.getOrbitalIrrep(l_in))
         {
            std::fstream fs;
            std::string filename = getenv("SAVE_H5_PATH");
            filename += "/orbs-scan-" + std::to_string(k_in) + "-" + std::to_string(l_in) + ".txt";
            fs.open(filename, std::fstream::out | std::fstream::trunc);

            fs.precision(10);

            fs << "# theta\trot\trot+v2dm" << std::endl;

            auto found = rdm.find_min_angle(k_in,l_in,0.3,getT,getV);

            std::cout << "Min:\t" << k_in << "\t" << l_in << "\t" << found.first << "\t" << found.second << std::endl;

            std::cout << "######################" << std::endl;


            int Na = 100;
#pragma omp parallel for
            for(int a=0;a<=Na;a++)
            {
               double theta = 1.0*M_PI/(1.0*Na) * a - M_PI/2.0;

               BoundaryPoint mymethod(ham);
               mymethod.set_tol_PD(1e-7);

//               mymethod.set_output(false);

               mymethod.getHam() = orig_ham;
               mymethod.getRDM() = rdm;
               mymethod.getHam().rotate(k_in, l_in, theta, getT, getV);

               double new_en = mymethod.evalEnergy();

               mymethod.Run();

#pragma omp critical
               fs << theta << "\t" << new_en << "\t" << mymethod.evalEnergy() << std::endl;
            }

            fs.close();
         }
}


/*  vim: set ts=3 sw=3 expandtab :*/
