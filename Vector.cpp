#include <iostream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <cstring>
#include <assert.h>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "Vector.h"
#include "lapack.h"

using namespace doci2DM;

/**
 * constructor that takes dimension as input
 * @param n dimension of the vector
 */
Vector::Vector(int n)
{
   this->n = n;

   vector.reset(new double[n]);
}

/**
 * Construct and initialize the Vector object by diagonalizing a Matrix object:
 */
Vector::Vector(Matrix &matrix)
{
   //allocate
   this->n = matrix.gn();

   vector.reset(new double[n]);

   diagonalize(matrix);
}

/**
 * copy constructor 
 * @param vec_copy The vector you want to be copied into the object you are constructing
 */
Vector::Vector(const Vector &vec_copy)
{
   this->n = vec_copy.n;

   vector.reset(new double[n]);

   std::memcpy(vector.get(), vec_copy.vector.get(), n*sizeof(double));
}

Vector::Vector(Vector &&vec_copy)
{
   this->n = vec_copy.n;

   vector = std::move(vec_copy.vector);
}


/**
 * Destructor
 */
Vector::~Vector()
{
}

/**
 * overload the equality operator
 * @param vector_copy The vector you want to be copied into this
 */
Vector &Vector::operator=(const Vector &vec_copy)
{
   assert(vec_copy.n == n);
   std::memcpy(vector.get(), vec_copy.vector.get(), n*sizeof(double));

   return *this;
}

Vector &Vector::operator=(Vector &&vec_copy)
{
   assert(vec_copy.n == n);

   vector = std::move(vec_copy.vector);

   return *this;
}

/**
 * Make all the number in your vector equal to the number a, e.g. usefull for initialization (Vector M = 0)
 * @param a the number
 */
Vector &Vector::operator=(double a)
{
   for(int i = 0;i < n;++i)
      vector[i] = a;

   return *this;
}

/**
 * overload the += operator for matrices
 * @param vector_pl The vector you want to add to this
 */
Vector &Vector::operator+=(const Vector &vector_pl)
{
   daxpy(1,vector_pl);

   return *this;
}

/**
 * overload the -= operator for matrices
 * @param vector_pl The vector you want to deduct from this
 */
Vector &Vector::operator-=(const Vector &vector_pl)
{
   daxpy(-1,vector_pl);

   return *this;
}

/**
 * add the vector vector_pl times the constant alpha to this
 * @param alpha the constant to multiply the vector_pl with
 * @param vector_pl the Vector to be multiplied by alpha and added to this
 */
Vector &Vector::daxpy(double alpha,const Vector &vector_pl)
{
   int inc = 1;

   daxpy_(&n,&alpha,vector_pl.vector.get(),&inc,vector.get(),&inc);

   return *this;
}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
Vector &Vector::operator/=(double c)
{
   dscal(1.0/c);

   return *this;
}

/**
 * *= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
Vector &Vector::operator*=(double c)
{
   dscal(c);

   return *this;
}

/**
 * write access to your vector, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &Vector::operator[](int i)
{
   assert(i<n);

   return vector[i];
}

/**
 * read access to your vector, change the number on index i: const version
 * @param i row number
 * @return the entry on place i
 */
double Vector::operator[](int i) const
{
   assert(i<n);

   return vector[i];
}

/**
 * Diagonalize the Matrix matrix when you have allready allocated the memory of the vector
 * on the correct dimension.
 */
void Vector::diagonalize(Matrix &matrix)
{
   *this = matrix.diagonalize();
}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications
 */
double *Vector::gVector()
{
   return vector.get();
}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications: const version
 */
const double *Vector::gVector() const
{
   return vector.get();
}

/**
 * @return the dimension of the vector
 */
int Vector::gn() const
{
   return n;
}

/**
 * @return the sum of all the elements in the vector
 */
double Vector::sum() const
{
   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += vector[i];

   return ward;
} 

double Vector::trace() const
{
   return sum();
}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double Vector::log_product() const
{
   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += log(vector[i]);

   return ward;
}

/**
 * @return inproduct of (*this) vector with vector_i
 * @param vector_i input vector
 */
double Vector::ddot(const Vector &vector_i) const
{
   int inc = 1;

   return ddot_(&n,vector.get(),&inc,vector_i.vector.get(),&inc);
}

/**
 * Scale the vector (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Vector::dscal(double alpha)
{
   int inc = 1;

   dscal_(&n,&alpha,vector.get(),&inc);
}

/**
 * Fill the vector with random numbers.
 */
void Vector::fill_Random()
{
   fill_Random(time(NULL));
}

/**
 * Fill the vector with random numbers.
 * @param seed the seed to use
 */
void Vector::fill_Random(int seed)
{
   srand(seed);

   for(int i = 0;i < n;++i)
      vector[i] = (double) rand()/RAND_MAX;
}

namespace doci2DM 
{
   std::ostream &operator<<(std::ostream &output,const doci2DM::Vector &vector_p)
   {
      for(int i = 0;i < vector_p.gn();++i)
         output << i << "\t" << vector_p[i] << endl;

      return output;
   }
}

/**
 * @return the minimal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::min() const
{
   double min = vector[0];

   for(int i=1;i<n;i++)
      if(vector[i]<min)
         min = vector[i];

   return min;
}

/**
 * @return the maximal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::max() const
{
   double max = vector[n-1];

   for(int i=n-2;i>=0;i--)
      if(vector[i]>max)
         max = vector[i];

   return max;
}

void Vector::invert()
{
   for(int i = 0;i < n;++i)
      vector[i] = 1.0/vector[i];
}

void Vector::sqrt(int option)
{
   if(option == 1)
      for(int i = 0;i < n;++i)
         vector[i] = std::sqrt(vector[i]);
   else 
      for(int i = 0;i < n;++i)
         vector[i] = 1.0/std::sqrt(vector[i]);
}

void Vector::symmetrize()
{
}

void Vector::L_map(const Vector &a,const Vector &b)
{
   assert(a.n == n && b.n == n);

   for(int i = 0;i < n;++i)
      vector[i] = a[i] * b[i] * a[i];
}

Vector &Vector::mprod(const Vector &x,const Vector &y)
{
   assert(x.n == n && y.n == n);

   for(int i = 0;i < n;++i)
      vector[i] = x[i] * y[i];

   return *this;
}

/**
 * Sort the vector, small to large
 */
void Vector::sort()
{
   std::sort(vector.get(), &vector.get()[n] );
}

void Vector::sep_pm(Vector &pos, Vector &neg)
{
   assert(pos.n == n && neg.n == n);
   pos = 0;
   neg = 0;

   for(int i=0;i<n;i++)
      if(vector[i] < 0)
         neg[i] = vector[i];
      else
         pos[i] = vector[i];
}

/* vim: set ts=3 sw=3 expandtab :*/
