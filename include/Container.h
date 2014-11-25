#ifndef CONTAINER_H
#define CONTAINER_H

#include <iostream>
#include <cstdlib>
#include <vector>

#include "BlockStructure.h"

namespace doci2DM
{

/**
 * Class to store combination of Blocks and Vectors
 */
class Container
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << blockmatrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << blockmatrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param blockmatrix_p the BlockMatrix you want to print
    */
   friend std::ostream &operator<<(std::ostream &output,const doci2DM::Container &container);

   public:

      //constructor
      Container(int n, int m);

      //copy constructor
      Container(const Container &);

      Container(Container &&);

      //destructor
      virtual ~Container() = default;

      void setMatrixDim(int,int,int);

      void setVectorDim(int,int,int);

      //overload equality operator
      Container &operator=(const Container &);

      Container &operator=(Container &&);

      Container &operator=(double);

      //overload += operator
      Container &operator+=(const Container &);

      //overload -= operator
      Container &operator-=(const Container &);

      Container &daxpy(double alpha,const Container &);

      Container &operator*=(double);

      Container &operator/=(double);

      Container &mprod(const Container &,const Container &);

      //easy to change the numbers
      double &operator()(int block,int i,int j);

      //easy to access the numbers
      double operator()(int block,int i,int j) const;

      double &operator()(int block,int i);

      //easy to access the numbers
      double operator()(int block,int i) const;

      Matrix & getMatrix(int);
      const Matrix & getMatrix(int) const;

      Vector & getVector(int);
      const Vector & getVector(int) const;

      BlockMatrix & getMatrices();
      const BlockMatrix & getMatrices() const;

      BlockVector & getVectors();
      const BlockVector & getVectors() const;

      int gnr() const;

      int gnMatrix() const;

      int gnVector() const;

      int gdimMatrix(int) const;

      int gdimVector(int) const;

      int gdegMatrix(int) const;

      int gdegVector(int) const;

      double trace() const;

      double ddot(const Container &) const;

      void invert();

      void dscal(double);

      void fill_Random();

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int);

      void L_map(const Container &,const Container &);

      void symmetrize();

      void sep_pm(Container &,Container &);

   private:

      std::unique_ptr<BlockMatrix> matrix;
      std::unique_ptr<BlockVector> vector;
};

}

#endif

/*  vim: set ts=3 sw=3 expandtab :*/
