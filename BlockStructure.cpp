#include <assert.h>

#include "Matrix.h"
#include "Vector.h"
#include "BlockStructure.h"

using namespace doci2DM;

/**
 * constructor: Watch out, the matrices themself haven't been allocated yet. 
 * Only the Blockmatrix itself is allocated and the array containing the dimensions, but not initialized.
 * @param blockmatrix.size() nr of blocks in the blockmatrix
 */
   template<class BlockType>
BlockStructure<BlockType>::BlockStructure(int nr)
{
   blocks.resize(nr);

   degen.resize(nr,0);
}

/**
 * copy constructor, make sure the input matrix and all the blocks have been allocated and filled before the copying
 * @param blockmat_copy The blockmatrix you want to be copied into the object you are constructing
 */
   template<class BlockType>
BlockStructure<BlockType>::BlockStructure(const BlockStructure<BlockType> &blockmat_copy)
{
   blocks.resize(blockmat_copy.blocks.size());

#pragma omp parallel for
   for(int i=0;i<blocks.size();i++)
      blocks[i].reset(new BlockType(blockmat_copy[i]));

   degen = blockmat_copy.degen;
}

   template<class BlockType>
BlockStructure<BlockType>::BlockStructure(BlockStructure<BlockType> &&blockmat_copy)
{
   blocks = std::move(blockmat_copy.blocks);

   degen = std::move(blockmat_copy.degen);
}

/**
 * function that allocates the memory of the blocks block with dimension dim sets and the degeneracy of the matrix
 * @param block the index of the block that will be allocated
 * @param dim the dimension of the particular block to be allocated
 * @param degeneracy the degeneracy of block "block"
 */
   template<class BlockType>
void BlockStructure<BlockType>::setDim(int block,int dim,int degeneracy)
{
   assert(block < blocks.size());

   this->degen[block] = degeneracy;

   blocks[block].reset(new BlockType(dim));
}

/**
 * [] overloaded, will return a reference to the block block in the blocks.
 * @param block index of the block to be returned
 * @return A reference to the Matrix object located on blocks[block].
 */
   template<class BlockType>
BlockType &BlockStructure<BlockType>::operator[](int block)
{
   return *blocks[block];
}

/**
 * [] overloaded, const version, will return a const reference to the block block in the blocks.
 * @param block index of the block to be returned
 * @return A reference to the Matrix object located on blocks[block].
 */
template<class BlockType>
const BlockType &BlockStructure<BlockType>::operator[](int block) const
{
   return *blocks[block];
}


/**
 * overload the equality operator: Make sure the blocks in both matrices have been allocated to the same dimensions and have the same degeneracy!
 * @param blockmat_copy The matrix you want to be copied into this
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator=(const BlockStructure<BlockType> &blockmat_copy)
{
   blocks.resize(blockmat_copy.blocks.size());

#pragma omp parallel for
   for(int i=0;i<blocks.size();i++)
      blocks[i].reset(new BlockType(blockmat_copy[i]));

   return *this;
}

/**
 * Make all the numbers in your blocks equal to the number a, e.g. usefull for initialization (BlockStructure M = 0)
 * @param a the number
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator=(double a)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      *blocks[i] = a;

   return *this;
}

/**
 * overload the += operator for matrices
 * @param blockmat_pl The matrix you want to add to this
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator+=(const BlockStructure<BlockType> &blockmat_pl)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      *blocks[i] += blockmat_pl[i];

   return *this;
}

/**
 * overload the -= operator for matrices
 * @param blockmat_pl The matrix you want to deduct from this
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator-=(const BlockStructure<BlockType> &blockmat_pl)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      *blocks[i] -= blockmat_pl[i];

   return *this;
}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param blockmat_pl the BlockStructure to be multiplied by alpha and added to this
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::daxpy(double alpha,const BlockStructure<BlockType> &blockmat_pl)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->daxpy(alpha,blockmat_pl[i]);

   return *this;
}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator/=(double c)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      *blocks[i] /= c;

   return *this;
}

   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::operator*=(double c)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      *blocks[i] *= c;

   return *this;
}

namespace doci2DM {
/**
 * write access to your blocks, change the number in block "block", on row i and column j
 * @param block The index of the block you want to access
 * @param i row number
 * @param j column number
 * @return the entry on place block,i,j
 */
   template<>
double &BlockStructure<Matrix>::operator()(int block,int i,int j)
{
   return (*blocks[block])(i,j);
}

/**
 * read access to your blocks, read the number in block "block" on row i and column j
 * @param block The index of the block you want to access
 * @param i row number
 * @param j column number
 * @return the entry on place block,i,j
 */
template<>
double BlockStructure<Matrix>::operator()(int block,int i,int j) const
{
   return (*blocks[block])(i,j);
}

   template<>
double& BlockStructure<Vector>::operator()(int block,int index)
{
   return (*blocks[block])[index];
}

template<>
double BlockStructure<Vector>::operator()(int block,int index) const
{
   return (*blocks[block])[index];
}

}

/**
 * @return the blocks.size() of blocks
 */
template<class BlockType>
int BlockStructure<BlockType>::gnr() const
{
   return blocks.size();
}

/**
 * @return the dimension of the Matrix on the block with index i
 */
template<class BlockType>
int BlockStructure<BlockType>::gdim(int i) const
{
   return blocks[i]->gn();
}

/**
 * @return the degeneracy of the block with index i
 */
template<class BlockType>
int BlockStructure<BlockType>::gdeg(int i) const
{
   return degen[i];
}

/**
 * @return the trace of the matrix, each block matrix is weighed with its degeneracy.
 */
template<class BlockType>
double BlockStructure<BlockType>::trace() const
{
   double ward = 0;

#pragma omp parallel for reduction(+:ward)
   for(int i = 0;i < blocks.size();++i)
      ward += degen[i]*blocks[i]->trace();

   return ward;
}

/**
 * @return inproduct of (*this) blocks with blocks_in, defined as Tr (A B)
 * @param blocks_in input matrix
 */
template<class BlockType>
double BlockStructure<BlockType>::ddot(const BlockStructure<BlockType> &blocks_in) const
{
   double ward = 0.0;

#pragma omp parallel for reduction(+:ward)
   for(int i = 0;i < blocks.size();i++)
      ward += degen[i]*blocks[i]->ddot(blocks_in[i]);

   return ward;
}

/**
 * Invert positive semidefinite symmetric blocks which is stored in (*this), original matrix (*this) is destroyed
 */
   template<class BlockType>
void BlockStructure<BlockType>::invert()
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->invert();
}

/**
 * Scale the blocks (*this) with parameter alpha
 * @param alpha scalefactor
 */
   template<class BlockType>
void BlockStructure<BlockType>::dscal(double alpha)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->dscal(alpha);
}

/**
 * Fill the matrix with random numbers.
 */
   template<class BlockType>
void BlockStructure<BlockType>::fill_Random()
{
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->fill_Random();
}

/**
 * Take the square root out of the positive semidefinite blocks, destroys original blocks, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
   template<class BlockType>
void BlockStructure<BlockType>::sqrt(int option)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->sqrt(option);
}

/**
 * Multiply symmetric blocks object left en right with symmetric blocks map to 
 * form another symmetric blocks and put it in (*this): this = map*object*map
 * @param map BlockStructure that will be multiplied to the left en to the right of matrix object
 * @param object central BlockStructure
 */
   template<class BlockType>
void BlockStructure<BlockType>::L_map(const BlockStructure<BlockType> &map,const BlockStructure<BlockType> &object)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->L_map(map[i],object[i]);
}

/**
 * BlockStructure product of two general blockmatrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
   template<class BlockType>
BlockStructure<BlockType> &BlockStructure<BlockType>::mprod(const BlockStructure<BlockType> &A, const BlockStructure<BlockType> &B)
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->mprod(A[i],B[i]);

   return *this;
}

/**
 * Copy upper triangle into lower triangle.
 */
   template<class BlockType>
void BlockStructure<BlockType>::symmetrize()
{
#pragma omp parallel for
   for(int i = 0;i < blocks.size();++i)
      blocks[i]->symmetrize();
}

namespace doci2DM {
template<>
void BlockStructure<Vector>::sort()
{
   for(unsigned int i = 0;i < blocks.size();++i)
      blocks[i]->sort();
}
}

   template<class BlockType>
void BlockStructure<BlockType>::sep_pm(BlockStructure<BlockType> &pos, BlockStructure<BlockType> &neg)
{
   for(unsigned int i=0;i<blocks.size();++i)
      blocks[i]->sep_pm(*pos.blocks[i], *neg.blocks[i]);
}


namespace doci2DM
{
   template<class MyBlockType>
   std::ostream &operator<< (std::ostream &output,const doci2DM::BlockStructure<MyBlockType> &blocks_p)
   {
      for(int i = 0;i < blocks_p.blocks.size();++i)
      {
         output << i << "\t" << blocks_p.blocks[i]->gn() << "\t" << blocks_p.degen[i] << std::endl;
         output << std::endl;

         output << *blocks_p.blocks[i] << std::endl;
      }

      return output;
   }

   // explicit instantiation must occur in namespace
   template class BlockStructure<doci2DM::Matrix>;
   template class BlockStructure<doci2DM::Vector>;

   // also instantiation for friend functions
   template std::ostream &operator<< (std::ostream &output,const doci2DM::BlockStructure<doci2DM::Matrix> &blocks_p);
   template std::ostream &operator<< (std::ostream &output,const doci2DM::BlockStructure<doci2DM::Vector> &blocks_p);
}


/* vim: set ts=3 sw=3 expandtab :*/
