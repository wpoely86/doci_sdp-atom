#include <assert.h>

#include "Container.h"

/**
 * @param n number of matrix
 * @param m number of vectors
 */
Container::Container(int n, int m)
{
   std::cout << "Fresh new container" << std::endl;
   matrix.reset(new BlockMatrix(n));

   vector.reset(new BlockVector(n));
}

Container::Container(const Container &orig)
{
   std::cout << "Copy Container" << std::endl;
   matrix.reset(new BlockMatrix(*orig.matrix));
   vector.reset(new BlockVector(*orig.vector));
}

Container::Container(Container &&orig)
{
   std::cout << "Moving Container" << std::endl;
   matrix = std::move(orig.matrix);
   vector = std::move(orig.vector);
}

Container::~Container()
{
   std::cout << "Bye Container" << std::endl;
}

/**
 * Set dimension and degeneracy for a Matrix
 * @param block which matrix
 * @param dim the dimensions
 * @param degen the degeneracy of the matrix
 */
void Container::setMatrixDim(int block,int dim,int degen)
{
   matrix->setDim(block, dim, degen);
}

/**
 * Set dimension and degeneracy for a Vector
 * @param block which vector
 * @param dim the dimensions
 * @param degen the degeneracy of the vector
 */
void Container::setVectorDim(int block,int dim,int degen)
{
   vector->setDim(block, dim, degen);
}

Container &Container::operator=(const Container &orig)
{
   std::cout << "Copy Assign Container" << std::endl;
   (*matrix) = *orig.matrix;
   (*vector) = *orig.vector;

   return *this;
}

Container &Container::operator=(Container &&orig)
{
   std::cout << "Move Assign Container" << std::endl;
   matrix = std::move(orig.matrix);
   vector = std::move(orig.vector);

   return *this;
}

Container &Container::operator=(double a)
{
   (*matrix) = a;
   (*vector) = a;

   return *this;
}

Container &Container::operator+=(const Container &orig)
{
   (*matrix) += *orig.matrix;
   (*vector) += *orig.vector;

   return *this;

}

Container &Container::operator-=(const Container &orig)
{
   (*matrix) -= *orig.matrix;
   (*vector) -= *orig.vector;

   return *this;
}

Container &Container::daxpy(double alpha,const Container &orig)
{
   matrix->daxpy(alpha, *orig.matrix);
   vector->daxpy(alpha, *orig.vector);

   return *this;
}

Container &Container::operator*=(double a)
{
   (*matrix) *= a;
   (*vector) *= a;

   return *this;
}

Container &Container::operator/=(double a)
{
   (*matrix) /= a;
   (*vector) /= a;

   return *this;
}

/**
 * General matrix matrix product of A and B. Stores result in *this.
 * @param A matrix A
 * @param B matrix B
 * @returns A*B (in *this)
 */
Container &Container::mprod(const Container &A,const Container &B)
{
   matrix->mprod(*A.matrix, *B.matrix);
   vector->mprod(*A.vector, *B.vector);

   return *this;
}

double &Container::operator()(int block,int i,int j)
{
   return (*matrix)(block,i,j);
}

double Container::operator()(int block,int i,int j) const
{
   return (*matrix)(block,i,j);
}

double &Container::operator()(int block,int i)
{
   return (*vector)(block,i);
}

double Container::operator()(int block,int i) const
{
   return (*vector)(block,i);
}

Matrix & Container::getMatrix(int i)
{
   return (*matrix)[i];
}

const Matrix & Container::getMatrix(int i) const
{
   return (*matrix)[i];
}

Vector & Container::getVector(int i)
{
   return (*vector)[i];
}

const Vector & Container::getVector(int i) const
{
   return (*vector)[i];
}

int Container::gnr() const
{
   return matrix->gnr() + vector->gnr();
}

int Container::gnMatrix() const
{
   return matrix->gnr();
}

int Container::gnVector() const
{
   return vector->gnr();
}

int Container::gdimMatrix(int i) const
{
   return matrix->gdim(i);
}

int Container::gdimVector(int i) const
{
   return vector->gdim(i);
}

int Container::gdegMatrix(int i) const
{
   return matrix->gdeg(i);
}

int Container::gdegVector(int i) const
{
   return vector->gdeg(i);
}

double Container::trace() const
{
   return matrix->trace() + vector->trace();
}

double Container::ddot(const Container &A) const
{
   return matrix->ddot(*A.matrix) + vector->ddot(*A.vector);
}

void Container::invert()
{
   matrix->invert();
   vector->invert();
}

void Container::dscal(double a)
{
   matrix->dscal(a);
   vector->dscal(a);
}

void Container::fill_Random()
{
   matrix->fill_Random();
   vector->fill_Random();
}

void Container::sqrt(int option)
{
   matrix->sqrt(option);
   vector->sqrt(option);
}

void Container::L_map(const Container &A,const Container &B)
{
   matrix->L_map(*A.matrix,*B.matrix);
   vector->L_map(*A.vector,*B.vector);
}

void Container::symmetrize()
{
   matrix->symmetrize();
}

/*  vim: set ts=3 sw=3 expandtab :*/
