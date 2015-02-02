#include <assert.h>

#include "include.h"

using namespace doci2DM;

// default empty
std::unique_ptr<helpers::tmatrix<int>> PHM::s2ph = nullptr;
std::unique_ptr<helpers::tmatrix<int>> PHM::ph2s = nullptr;
std::unique_ptr<helpers::tmatrix<int>> PHM::s2b = nullptr;
std::unique_ptr<helpers::tmatrix<int>> PHM::b2s = nullptr;

/**
 * @param L the number of levels
 * @param N the number of particles
 */
PHM::PHM(int L, int N): BlockMatrix(1+(L*(L-1))/2)
{
   this->L = L;
   this->N = N;

   if(!s2ph || !ph2s)
      constr_lists(L);

   // one LxL block
   setDim(0, L, 1);

   // all the rest are 2x2 blocks
   for(int i=1;i<gnr();i++)
      setDim(i, 2, 1);
}

void PHM::constr_lists(int L)
{
   int M = 2*L;
   int n_ph = M*M;

   s2ph.reset(new helpers::tmatrix<int>(M,M));
   (*s2ph) = -1; // if you use something you shouldn't, this will case havoc

   ph2s.reset(new helpers::tmatrix<int>(n_ph,2));
   (*ph2s) = -1; // if you use something you shouldn't, this will case havoc

   int tel = 0;

   // a b
   for(int a=0;a<L;a++)
      for(int b=0;b<L;b++)
         (*s2ph)(a,b) = tel++;

   // \bar a \bar b
   for(int a=L;a<M;a++)
      for(int b=L;b<M;b++)
         (*s2ph)(a,b) = tel++;

   // a \bar b 
   for(int a=0;a<L;a++)
      for(int b=L;b<M;b++)
         (*s2ph)(a,b) = tel++;

   // \bar a b
   for(int a=L;a<M;a++)
      for(int b=0;b<L;b++)
         (*s2ph)(a,b) = tel++;

   assert(tel == n_ph);

   for(int a=0;a<M;a++)
      for(int b=0;b<M;b++)
      {
         (*ph2s)((*s2ph)(a,b),0) = a;
         (*ph2s)((*s2ph)(a,b),1) = b;
      }

   // these convert between the 2x2 block index and the sp index for that block
   s2b.reset(new helpers::tmatrix<int>(L,L));
   (*s2b) = -1; // if you use something you shouldn't, this will case havoc

   b2s.reset(new helpers::tmatrix<int>((L*(L-1))/2+1,2));
   (*b2s) = -1; // if you use something you shouldn't, this will case havoc

   // sp index to block index
   tel = 1;
   for(int a=0;a<L;a++)
      for(int b=a+1;b<L;b++)
      {
         (*s2b)(a,b) = (*s2b)(b,a) = tel;
         (*b2s)(tel,0) = a;
         (*b2s)(tel,1) = b;
         ++tel;
      }

   assert(tel == b2s->getn());
}

/**
 * The access operator for the G matrix. It's a bit complicated due
 * to the fact that we don't store the full blown G matrix. We construct
 * each element on the fly
 * @param a the first sp index
 * @param b the second sp index
 * @param c the thirth sp index
 * @param d the fourth sp index
 * @return the G matrix element
 */
double PHM::operator()(int a, int b, int c, int d) const
{
   double res = 0;

   auto getelem = [this](int a, int b) {
      if(a!=b)
      {
         int i = (*s2b)(a,b);
         int j = a > b ? 1 : 0;
         return (*this)(i,j,j);
      } else
         return (*this)(0,a,a);
   };

   auto getrho = [this](int a, int b) {
      if(a!=b)
      {
         int i = (*s2b)(a,b);
         return -1 * (*this)(i,0,1);
      } else
         return (*this)(0,a,a);
   };

   // \bar a b ; \bar c d 
   // a \bar b ; c \bar d 
   if(a/L != b/L && c/L != d/L && a/L==c/L)
   {
      if(a%L==c%L && b%L==d%L)
         res += getelem(a%L,b%L);

      if(a%L==d%L && b%L==c%L)
         res -= getrho(a%L,b%L);
   }
   // check that a and b are both up or down (same for c and d)
   else if(a/L == b/L && c/L == d/L)
   {
      // if a and c have the same spin
      if(a/L == c/L)
      {
         if(a==c && b==d)
            res += getelem(a%L,b%L);

         if(a==b && c==d && a!=c)
            res += (*this)(0,a%L,c%L);
      } else
      {
         if(a%L==d%L && b%L==c%L)
            res += getrho(a%L,b%L);

         if(a==b && c==d && a%L!=c%L)
            res += (*this)(0,a%L,c%L);
      }
   }

   return res;
}

int PHM::gN() const
{
   return N;
}

int PHM::gL() const
{
   return L;
}

/**
 * Get a 2x2 block of the G matrix, corresponding with
 * indices a and b
 * @param a the first sp index
 * @param b the second sp index
 * @return the 2x2 matrix
 */
const Matrix& PHM::getBlock(int a, int b) const
{
   const int idx = (*s2b)(a,b);
   assert(idx>0);

   return (*this)[idx];
}

/**
 * Get a 2x2 block of the G matrix, corresponding with
 * indices a and b
 * @param a the first sp index
 * @param b the second sp index
 * @return the 2x2 matrix
 */
Matrix& PHM::getBlock(int a, int b)
{
   const int idx = (*s2b)(a,b);
   assert(idx>0);

   return (*this)[idx];
}

/**
 * Create the G condition from a TPM
 * @param tpm the TPM to use
 */
void PHM::G(const TPM &tpm)
{
   const SPM spm(tpm);

   // first the LxL block
   for(int a=0;a<L;a++)
   {
      for(int b=a;b<L;b++)
         (*this)(0,b,a) = (*this)(0,a,b) = tpm.getDiag(a,b);

      (*this)(0,a,a) += spm(0,a);
   }

   // now all the 2x2 blocks
   for(int i=1;i<gnr();i++)
      {
         auto& block = (*this)[i];
         int a = (*b2s)(i,0);
         int b = (*b2s)(i,1);

         block(0,0) = spm(0,a) - tpm.getDiag(a,b);
         block(1,1) = spm(0,b) - tpm.getDiag(a,b);
         block(0,1) = block(1,0) = - tpm(a,a+L,b,b+L);
      }
}

/**
 * Build the G matrix using the full
 * G matrix image
 * @param tpm the TPM to use
 * @return the full G matrix
 */
Matrix PHM::Gimg(const TPM &tpm) const
{
   int M = 2*L;
   Matrix Gmat(M*M);

   SPM spm(tpm);

   for(int i=0;i<M*M;++i)
   {
      int a = (*ph2s)(i,0);
      int b = (*ph2s)(i,1);

      for(int j=i;j<M*M;++j)
      {
         int c = (*ph2s)(j,0);
         int d = (*ph2s)(j,1);

         Gmat(i,j) = -tpm(a,d,c,b);

         if(b == d && a == c)
            Gmat(i,j) += spm(0,a%L);

         Gmat(j,i) = Gmat(i,j);
      }
   }

   return Gmat;
}

/**
 * Build the full G matrix out of *this
 * object
 * @return the full dense G matrix
 */
Matrix PHM::Gbuild() const
{
   const int M = 2*L;
   Matrix Gmat(M*M);

   for(int i=0;i<M*M;++i)
   {
      int a = (*ph2s)(i,0);
      int b = (*ph2s)(i,1);

      for(int j=i;j<M*M;++j)
      {
         int c = (*ph2s)(j,0);
         int d = (*ph2s)(j,1);

         Gmat(i,j) = Gmat(j,i) = (*this)(a,b,c,d);
      }
   }

   return Gmat;
}

namespace doci2DM 
{
   std::ostream &operator<<(std::ostream &output,doci2DM::PHM &phm)
   {
      output << "The LxL block: " << std::endl;
      output << phm[0] << std::endl;

      output << std::endl;

      for(int i=1;i<phm.gnr();i++)
      {
         output << "Block " << i-1 << " for " << (*phm.b2s)(i,0) << "\t" << (*phm.b2s)(i,1) << std::endl;
         output << phm[i] << std::endl;
      }

      return output;
   }
}

/**
 * Use the 2x2 sep_pm for all the 2x2 blocks instead
 * of the full blown sep_pm
 * @param pos the positive part of the matrix
 * @param neg the negative part of the matrix
 */
void PHM::sep_pm(BlockMatrix &pos, BlockMatrix &neg)
{
   (*this)[0].sep_pm(pos[0], neg[0]);

#pragma omp parallel for if(gnr()>100)
   for(int i=1;i<gnr();i++)
      (*this)[i].sep_pm_2x2(pos[i], neg[i]);
}

/**
 * Pull the sqrt out of this matrix.
 * Optimized for the structure of the 2x2 matrices
 * @param option 1 => positive sqrt, -1 => negative sqrt
 */
void PHM::sqrt(int option)
{
   (*this)[0].sqrt(option);

#pragma omp parallel for if(gnr()>100)
   for(int i=1;i<gnr();i++)
      (*this)[i].sqrt_2x2(option);
}

void PHM::invert()
{
   (*this)[0].invert();

#pragma omp parallel for if(gnr()>100)
   for(int i=1;i<gnr();i++)
      (*this)[i].invert_2x2();
}

void PHM::L_map(const BlockMatrix &map,const BlockMatrix &object)
{
   (*this)[0].L_map(map[0], object[0]);

#pragma omp parallel for if(gnr()>100)
   for(int i=1;i<gnr();i++)
      (*this)[i].L_map_2x2(map[i], object[i]);
}

void PHM::WriteToFile(hid_t &group_id) const
{
   hid_t       dataset_id, dataspace_id;
   herr_t      status;

   // first the LxL block
   hsize_t dimblock = (*this)[0].gn() * (*this)[0].gn();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Block", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   double *data = const_cast<Matrix &>((*this)[0]).gMatrix();

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   // first the 2x2 block
   dimblock = 4;
   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   for(int i=1;i<gnr();i++)
   {
      std::string blockname = "2x2_" + std::to_string(i);

      dataset_id = H5Dcreate(group_id, blockname.c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      data = const_cast<Matrix &>((*this)[i]).gMatrix();

      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      HDF5_STATUS_CHECK(status);

      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
   }

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);
}


void PHM::WriteFullToFile(hid_t &group_id) const
{
   hid_t       dataset_id, attribute_id, dataspace_id;
   herr_t      status;

   int M = 2*L;

   Matrix fullTPM(M*M);
   fullTPM = 0;

   for(int a=0;a<M;a++)
      for(int b=0;b<M;b++)
         for(int c=0;c<M;c++)
            for(int d=0;d<M;d++)
            {
               int idx1 = (*s2ph)(a,b);
               int idx2 = (*s2ph)(c,d);

               if(idx1>=0 && idx2>=0)
                  fullTPM(idx1, idx2) = (*this)(a,b,c,d);
            }


   hsize_t dimblock = M*M*M*M;

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "PHM", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fullTPM.gMatrix());
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   dataspace_id = H5Screate(H5S_SCALAR);

   attribute_id = H5Acreate (dataset_id, "M", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &M );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Acreate (dataset_id, "N", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &N );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);
}

/*  vim: set ts=3 sw=3 expandtab :*/
