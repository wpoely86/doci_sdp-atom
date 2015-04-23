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
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @END LICENSE
 */

#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <assert.h>
#include <hdf5.h>

using std::endl;
using std::ostream;

#include "Matrix.h"
#include "lapack.h"
#include "Vector.h"

#define HDF5_STATUS_CHECK(status) {                 \
    if(status < 0)                                  \
    std::cerr << __FILE__ << ":" << __LINE__ <<     \
    ": Problem with writing to file. Status code="  \
    << status << std::endl;                         \
}


using namespace doci2DM;

/**
 * constructor 
 * @param n dimension of the matrix
 */
Matrix::Matrix(int n)
{
   this->n = n;

   // we store column major (for Fortran compatiblity)
   matrix.reset(new double[n*n]);
}

/**
 * copy constructor 
 * @param mat_copy The matrix you want to be copied into the object you are constructing
 */
Matrix::Matrix(const Matrix &orig){

   this->n = orig.n;

   matrix.reset(new double[n*n]);

   std::memcpy(matrix.get(), orig.matrix.get(), n*n*sizeof(double));
}

Matrix::Matrix(Matrix &&orig)
{
   n = orig.n;
   matrix = std::move(orig.matrix);
}

/**
 * Destructor
 */
Matrix::~Matrix()
{
}

/**
 * overload the equality operator
 * @param matrix_copy The matrix you want to be copied into this
 */
Matrix &Matrix::operator=(const Matrix &matrix_copy)
{
   assert(n == matrix_copy.n);

   std::memcpy(matrix.get(), matrix_copy.matrix.get(), n*n*sizeof(double));

   return *this;
}

Matrix &Matrix::operator=(Matrix &&matrix_copy)
{
   assert(n == matrix_copy.n);

   matrix = std::move(matrix_copy.matrix);

   return *this;
}

/**
 * Make all the number in your matrix equal to the number a, e.g. usefull for initialization (Matrix M = 0)
 * @param a the number
 */
Matrix &Matrix::operator=(double a){

   for(int i = 0;i < n*n;++i)
      matrix[i] = a;

   return *this;
}

/**
 * overload the += operator for matrices
 * @param matrix_pl The matrix you want to add to this
 */
Matrix &Matrix::operator+=(const Matrix &matrix_pl)
{
   int dim = n*n;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix.get(),&inc,matrix.get(),&inc);

   return *this;
}

/**
 * overload the -= operator for matrices
 * @param matrix_pl The matrix you want to deduct from this
 */
Matrix &Matrix::operator-=(const Matrix &matrix_pl)
{
   int dim = n*n;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix.get(),&inc,matrix.get(),&inc);

   return *this;
}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param matrix_pl the Matrix to be multiplied by alpha and added to this
 */
Matrix &Matrix::daxpy(double alpha,const Matrix &matrix_pl)
{
   int dim = n*n;
   int inc = 1;

   daxpy_(&dim,&alpha,matrix_pl.matrix.get(),&inc,matrix.get(),&inc);

   return *this;
}

/**
 * *= operator overloaded: multiply by a constant
 * @param c the number to multiply your matrix with
 */
Matrix &Matrix::operator*=(double c)
{
   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&c,matrix.get(),&inc);

   return *this;
}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
Matrix &Matrix::operator/=(double c)
{
   operator*=(1.0/c);

   return *this;
}

/**
 * write access to your matrix, change the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double &Matrix::operator()(int i,int j)
{
   assert(i<n && j<n);

   return matrix[i+j*n];
}

/**
 * read access to your matrix, view the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double Matrix::operator()(int i,int j) const
{
   assert(i<n && j<n);

   return matrix[i+j*n];
}

/**
 * @return the underlying pointer to matrix, useful for mkl applications
 */
double *Matrix::gMatrix()
{
   return matrix.get();
}

const double *Matrix::gMatrix() const
{
   return matrix.get();
}

/**
 * @return the dimension of the matrix
 */
int Matrix::gn() const
{
   return n;
}

/**
 * @return the trace of the matrix:
 */
double Matrix::trace() const
{
   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += matrix[i+i*n];

   return ward;
}

/**
 * Diagonalizes symmetric matrices. Watch out! The current matrix (*this) is destroyed, in it
 * the eigenvectors will be stored (one in every column).
 * @param eigenvalues the pointer of doubles in which the eigenvalues will be storen, Watch out, its memory
 * has to be allocated on the dimension of the matrix before you call the function.
 */
Vector Matrix::diagonalize()
{
   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   Vector eigenvalues(n);

   std::unique_ptr<double []> work (new double [lwork]);

   int info = 0;

   dsyev_(&jobz,&uplo,&n,matrix.get(),&n,eigenvalues.gVector(),work.get(),&lwork,&info);

   if(info)
      std::cerr << "dsyev failed. info = " << info << std::endl;

   return eigenvalues;
}

/**
 * Diagonalize 2x2 matrix
 * Overwrite matrix with the (normed) eigenvectors
 * @return Vector with the eigenvalues
 */
Vector Matrix::diagonalize_2x2()
{
   assert(n==2 && "Only works for 2x2 matrices");

   // matrix has this form:
   // a c
   // c d
   auto a = matrix[0];
   auto c = matrix[1];
   auto &d = matrix[3];

   // double discr = (a+d)*(a+d) - 4*(a*d-c*c);
   double discr = a*a + d*d + 4*c*c - 2*a*d;
   // sometimes, we get -1e15. Deal with it.
   if(discr<0)
      discr = 0;
//   assert(!(discr < 0) && "Impossible!?!");

   Vector eig(2);
   eig[0] = ((a+d) - std::sqrt(discr))/2;
   eig[1] = ((a+d) + std::sqrt(discr))/2;

   double norm = 1.0/std::sqrt(1 + (eig[0]-a)*(eig[0]-a)/(c*c));

   matrix[0] = norm;
   matrix[1] = norm*(eig[0]-a)/c;

   norm = 1.0/std::sqrt(1 + (eig[1]-a)*(eig[1]-a)/(c*c));

   matrix[2] = norm;
   matrix[3] = norm*(eig[1]-a)/c;

   return eig;
}

/**
 * @return inproduct of (*this) matrix with matrix_i, defined as Tr (A B)
 * @param matrix_i input matrix
 */
double Matrix::ddot(const Matrix &matrix_i) const
{
   int dim = n*n;
   int inc = 1;

   return ddot_(&dim,matrix.get(),&inc,matrix_i.matrix.get(),&inc);
}

/**
 * Invert positive semidefinite symmetric matrix is stored in (*this), original matrix (*this) is destroyed
 */
void Matrix::invert()
{
   char uplo = 'U';

   int info = 0;

   dpotrf_(&uplo,&n,matrix.get(),&n,&info);//cholesky decompositie

   if(info)
      std::cerr << "dpotrf failed. info = " << info << std::endl;

   dpotri_(&uplo,&n,matrix.get(),&n,&info);//inverse berekenen

   if(info)
      std::cerr << "dpotri failed. info = " << info << std::endl;

   //terug symmetrisch maken:
   this->symmetrize();
}

/**
 * Invert 2x2 matrix
 */
void Matrix::invert_2x2()
{
   assert(n==2 && "Only works for 2x2 matrices");

   // matrix has this form:
   // a c
   // c d
   auto a = matrix[0];
   auto &c = matrix[1];
   auto &d = matrix[3];

   double fac = 1.0/(a*d-c*c);

   matrix[0] = fac * d;
   matrix[3] = fac * a;
   matrix[1] *= -fac;
   matrix[2] *= -fac;
}

/**
 * Scale the matrix (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Matrix::dscal(double alpha)
{
   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&alpha,matrix.get(),&inc);
}

/**
 * Fill the matrix with random numbers. The seed is the current time.
 */
void Matrix::fill_Random()
{
   fill_Random(time(NULL));
}

/**
 * Fill the matrix with random numbers.
 * @param seed the seed for the random number generator
 */
void Matrix::fill_Random(int seed)
{
   srand(seed);

   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         matrix[i*n+j] = matrix[j*n+i] = (double) rand()/RAND_MAX;
}

/**
 * Take the square root out of the positive semidefinite matrix, destroys original matrix, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void Matrix::sqrt(int option)
{
   Matrix hulp(*this);

   auto eigen = hulp.diagonalize();

   if(option == 1)
      for(int i = 0;i < n;++i)
         eigen[i] = std::sqrt(eigen[i]);
   else
      for(int i = 0;i < n;++i)
         eigen[i] = 1.0/std::sqrt(eigen[i]);

   //hulp opslaan
   Matrix hulp_c = hulp;

   //vermenigvuldigen met diagonaalmatrix
   hulp_c.mdiag(eigen);

   //en tenslotte de laatste matrixvermenigvuldiging
   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&transA,&transB,&n,&n,&n,&alpha,hulp_c.matrix.get(),&n,hulp.matrix.get(),&n,&beta,matrix.get(),&n);
}

/**
 * Take the square root out of the positive semidefinite matrix, destroys original matrix, square root will be put in (*this)
 * Only works for 2x2 matrices
 * @param option = 1, positive square root, = -1, negative square root.
 */
void Matrix::sqrt_2x2(int option)
{
   assert(n==2 && "Only works for 2x2 matrices");

   // matrix has this form:
   // a c
   // c d
   auto a = matrix[0];
   auto c = matrix[1];
   auto &d = matrix[3];

   // double discr = (a+d)*(a+d) - 4*(a*d-c*c);
   double discr = a*a + d*d + 4*c*c - 2*a*d;
   // sometimes, we get -1e15. Deal with it.
   if(discr<0)
      discr = 0;
//   assert(!(discr < 0) && "Impossible!?!");

   // the eigenvalues
   double x1 = ((a+d) - std::sqrt(discr))/2;
   double x2 = ((a+d) + std::sqrt(discr))/2;

   double norm1 = 1.0/(1 + (x1-a)*(x1-a)/(c*c));
   double norm2 = 1.0/(1 + (x2-a)*(x2-a)/(c*c));

   // calculate the sqrt from both eigenvalues
   double x1_s, x2_s;

   if(option)
   {
      x1_s = std::sqrt(x1);
      x2_s = std::sqrt(x2);
   } else {
      x1_s = 1.0/std::sqrt(x1);
      x2_s = 1.0/std::sqrt(x2);
   }

   matrix[0] = norm1*x1_s + norm2*x2_s;
   matrix[1] = matrix[2] = norm1*x1_s*(x1-a)/c + norm2*x2_s*(x2-a)/c;
   matrix[3] = norm1*x1_s*(x1-a)/c*(x1-a)/c + norm2*x2_s*(x2-a)/c*(x2-a)/c;
}

/**
 * Multiply this matrix with diagonal matrix
 * @param diag Diagonal matrix to multiply with this, has to be allocated on matrix dimension.
 */
void Matrix::mdiag(const Vector &diag)
{
   int inc = 1;

   for(int i = 0;i < n;++i)
      dscal_(&n, &diag.gVector()[i], &matrix[i*n], &inc);
}

/**
 * Multiply symmetric matrix object left en right with symmetric matrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map matrix that will be multiplied to the left en to the right of matrix object
 * @param object central matrix
 */
void Matrix::L_map(const Matrix &map,const Matrix &object)
{
   char side = 'L';
   char uplo = 'U';

   double alpha = 1.0;
   double beta = 0.0;

   Matrix hulp(n);

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix.get(),&n,object.matrix.get(),&n,&beta,hulp.matrix.get(),&n);

   side = 'R';

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix.get(),&n,hulp.matrix.get(),&n,&beta,matrix.get(),&n);

   //expliciet symmetriseren van de uit matrix
   this->symmetrize();
}

/**
 * L_map for 2x2 matrices. Will store map*object*map in this
 * @param map matrix that will be multiplied to the left en to the right of matrix object
 * @param object central matrix
 */
void Matrix::L_map_2x2(const Matrix &map,const Matrix &object)
{
   assert(n==2 && map.n==2 && object.n==2 && "Only works for 2x2 matrices");

   // matrix has this form:
   // a c
   // c d
   auto &a1 = object.matrix[0];
   auto &c1 = object.matrix[1];
   auto &d1 = object.matrix[3];

   auto &a2 = map.matrix[0];
   auto &c2 = map.matrix[1];
   auto &d2 = map.matrix[3];

   matrix[0] = a2*(a1*a2+c1*c2)+c2*(a2*c1+c2*d1);
   matrix[1] = matrix[2] = a2*(a1*c2+c1*d2)+c2*(c1*c2+d1*d2);
   matrix[3] = c2*(a1*c2+c1*d2)+d2*(c1*c2+d1*d2);
}

/**
 * Matrix product of two general matrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
Matrix &Matrix::mprod(const Matrix &A,const Matrix &B)
{
   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&trans,&trans,&n,&n,&n,&alpha,A.matrix.get(),&n,B.matrix.get(),&n,&beta,matrix.get(),&n);

   return *this;
}

/**
 * Copy upper triangle into lower triangle.
 */
void Matrix::symmetrize()
{
   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[j+i*n] = matrix[i+n*j];
}

namespace doci2DM {
   std::ostream &operator<<(std::ostream &output,const doci2DM::Matrix &matrix_p)
   {
      //   output << std::setprecision(2) << std::fixed;
      //
      //   for(int i = 0;i < matrix_p.gn();i++)
      //   {
      //      for(int j = 0;j < matrix_p.gn();j++)
      //         output << std::setfill('0') << std::setw(6) << matrix_p(i,j) << " ";
      //
      //      output << endl;
      //   }
      //
      //   output << endl;
      //   output << std::setprecision(10) << std::scientific;

      for(int i = 0;i < matrix_p.gn();++i)
         for(int j = 0;j < matrix_p.gn();++j)
            output << i << "\t" << j << "\t" << matrix_p(i,j) << endl;

      //   output.unsetf(std::ios_base::floatfield);

      return output;
   }
}

void Matrix::SaveRawToFile(const std::string filename) const
{
    hid_t       file_id, dataset_id, dataspace_id, attribute_id;
    herr_t      status;

    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_STATUS_CHECK(file_id);

    hsize_t dimarr = n*n;

    dataspace_id = H5Screate_simple(1, &dimarr, NULL);

    dataset_id = H5Dcreate(file_id, "matrix", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix.get() );
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    dataspace_id = H5Screate(H5S_SCALAR);

    attribute_id = H5Acreate (dataset_id, "n", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite (attribute_id, H5T_NATIVE_INT, &n );
    HDF5_STATUS_CHECK(status);
    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Fclose(file_id);
    HDF5_STATUS_CHECK(status);
}

/**
 * Seperate matrix into two matrices, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void Matrix::sep_pm(Matrix &p,Matrix &m)
{
   //init:
#ifdef __INTEL_COMPILER
   p = 0;
   m = *this;
#else
   p = *this;
   m = 0;
#endif

   std::unique_ptr<double []> eigenvalues(new double [n]);

   //diagonalize orignal matrix:
   char jobz = 'V';
   char uplo = 'U';

   // to avoid valgrind errors
   int info = 0;

#ifdef SYEVD
   int lwork = 2*n*n+6*n+1;
   int liwork = 5*n+3;

   std::unique_ptr<double []> work(new double [lwork]);
   std::unique_ptr<int []> iwork(new int [liwork]);

   dsyevd_(&jobz,&uplo,&n,matrix.get(),&n,eigenvalues.get(),work.get(),&lwork,iwork.get(),&liwork,&info);

   if(info)
      std::cerr << "dsyevd in sep_pm failed..." << std::endl;
#else
   int lwork = 3*n - 1;

   std::unique_ptr<double []> work(new double [lwork]);

   dsyev_(&jobz,&uplo,&n,matrix.get(),&n,eigenvalues.get(),work.get(),&lwork,&info);

   if(info)
      std::cerr << "dsyev in sep_pm failed..." << std::endl;
#endif

#ifdef __INTEL_COMPILER
   if( eigenvalues[0] >= 0 )
   {
      p = m;
      m = 0;
      return;
   }

   if( eigenvalues[n-1] < 0)
      return;

   int i = 0;

   while(i < n && eigenvalues[i] < 0.0)
      eigenvalues[i++] = 0;

   Matrix copy(*this);

   int inc = 1;

   // multiply each column with an eigenvalue
   for(i = 0;i < n;++i)
      dscal_(&n,&eigenvalues[i],&matrix[i*n],&inc);

   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&transA,&transB,&n,&n,&n,&alpha,matrix.get(),&n,copy.matrix.get(),&n,&beta,p.matrix.get(),&n);

   m -= p;
#else
   if( eigenvalues[0] >= 0 )
      return;

   if( eigenvalues[n-1] < 0)
   {
      m = p;
      p = 0;
      return;
   }

   //fill the plus and minus matrix
   int i = 0;

   while(i < n && eigenvalues[i] < 0.0)
   {
      for(int j = 0;j < n;++j)
         for(int k = j;k < n;++k)
            m(j,k) += eigenvalues[i] * (*this)(j,i) * (*this)(k,i);

      ++i;
   }

   m.symmetrize();

   p -= m;
#endif
}

/**
 * Special version of sep_pm that only works for 2x2 matrices
 * @param pos the positive part of the matrix
 * @param neg the negative part of the matrix
 */
void Matrix::sep_pm_2x2(Matrix &pos,Matrix &neg)
{
   assert(n==2 && "Only works for 2x2 matrices");

   // matrix has this form:
   // a c
   // c d
   auto &a = matrix[0];
   auto &c = matrix[1];
   auto &d = matrix[3];

   // there are 3 cases: both positive, both negative and one of each

   // double discr = (a+d)*(a+d) - 4*(a*d-c*c);
   double discr = a*a + d*d + 4*c*c - 2*a*d;
   // sometimes, we get -1e15. Deal with it.
   if(discr<0)
      discr = 0;
//   assert(!(discr < 0) && "Impossible!?!");

   double x1 = (a+d) - std::sqrt(discr);
   double x2 = (a+d) + std::sqrt(discr);

   // both positive
   if(x1>0)
   {
      pos = *this;
      neg = 0;
      return;
   }

   // both negative
   if(x2<0)
   {
      pos = 0;
      neg = *this;
      return;
   }

   // one negative eigenvalue
   x1 /= 2; // the actual eigenvalue still have to be divided by 2
   double norm = 1.0/(1 + (x1-a)*(x1-a)/(c*c));

   neg(0,0) = x1*norm;
   neg(1,0) = neg(0,1) = x1*norm*(x1-a)/c;
   neg(1,1) = x1*norm*(x1-a)*(x1-a)/(c*c);

   pos = *this;
   pos -= neg;
}

/**
 * Set this matrix to a unit matrix
 */
void Matrix::unit()
{
   (*this) = 0;

   for(int i=0;i<n;i++)
      matrix[i*n+i] = 1.0;
}

/* vim: set ts=3 sw=3 expandtab :*/
