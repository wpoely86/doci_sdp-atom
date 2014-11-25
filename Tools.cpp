#include <hdf5.h>
#include "include.h"

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

   group_id = H5Gopen(file_id, "/integrals", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   attribute_id = H5Aopen(group_id, "nelectrons", H5P_DEFAULT);
   HDF5_STATUS_CHECK(attribute_id);

   status = H5Aread(attribute_id, H5T_NATIVE_INT, &nelectrons);
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
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

/*  vim: set ts=3 sw=3 expandtab :*/
