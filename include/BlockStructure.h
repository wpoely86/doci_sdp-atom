#ifndef BLOCKSTRUCTURE_H
#define BLOCKSTRUCTURE_H

#include <iostream>
#include <memory>
#include <vector>

#include "Matrix.h"
#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 5-04-2010\n\n
 * This is a class written for symmetric block matrices. It containts an array of Matrix objects and an 
 * array containing the dimensions of the different block. It redefines all the member functions of the Matrix
 * class, which uses the lapack and blas routines for matrix computations.
 */
template<class BlockType>
class BlockStructure
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

   friend std::ostream &operator<< (std::ostream &output,const BlockStructure<BlockType> &blocks_p)
   {
      for(int i = 0;i < blocks_p.blocks.size();++i)
      {
         output << i << "\t" << blocks_p.blocks[i]->gn() << "\t" << blocks_p.degen[i] << std::endl;
         output << std::endl;

         output << *blocks_p.blocks[i] << std::endl;
      }

      return output;
   }

   public:

      BlockStructure(int);

      //construct with as input a BlockMatrixType
      //   BlockStructure(BlockType &);

      //copy constructor
      BlockStructure(const BlockStructure<BlockType> &);

      BlockStructure(BlockStructure<BlockType> &&);

      //destructor
      virtual ~BlockStructure() = default;

      void setDim(int,int,int);

      //overload equality operator
      BlockStructure &operator=(const BlockStructure<BlockType> &);

      BlockStructure &operator=(BlockStructure<BlockType> &&);

      BlockStructure &operator=(double);

      //overload += operator
      BlockStructure &operator+=(const BlockStructure<BlockType> &);

      //overload -= operator
      BlockStructure &operator-=(const BlockStructure<BlockType> &);

      BlockStructure &daxpy(double alpha,const BlockStructure<BlockType> &);

      BlockStructure &operator*=(double);

      BlockStructure &operator/=(double);

      BlockType &operator[](int block);

      //const version
      const BlockType &operator[](int block) const;

      double &operator()(int block,int i, int j);

      //easy to access the numbers
      double operator()(int block,int i, int j) const;

      //easy to change the numbers
      double &operator()(int block,int index);

      //easy to access the numbers
      double operator()(int block,int index) const;

      void diagonalize(BlockStructure<Matrix> &);

      int gnr() const;

      int gdim(int) const;

      int gdeg(int) const;

      double trace() const;

      double ddot(const BlockStructure<BlockType> &) const;

      void dscal(double alpha);

      double min() const;

      double max() const;

      void fill_Random();

      void invert();

      void sqrt(int);

      void L_map(const BlockStructure<BlockType> &map,const BlockStructure<BlockType> &object);

      BlockStructure<BlockType> &mprod(const BlockStructure<BlockType> &A, const BlockStructure<BlockType> &B);

      void symmetrize();

      void sort();

   private:

      std::vector< std::unique_ptr<BlockType> > blocks;

      //!degeneracy of the blocks
      std::vector<int> degen;
};

typedef BlockStructure<Matrix> BlockMatrix;
typedef BlockStructure<Vector> BlockVector;

#endif

/* vim: set ts=3 sw=3 expandtab :*/
