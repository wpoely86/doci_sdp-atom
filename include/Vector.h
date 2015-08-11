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
 * v2DM-DOCI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with v2DM-DOCI.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @END LICENSE
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>

#include "Matrix.h"


namespace doci2DM
{

/**
 * @author Brecht Verstichel
 * @date 15-04-2010\n\n
 * This is a class written for vectors. It will contain the eigenvalues of the TPM, etc. Matrices. It is a template class,
 * corresponding to the different VectorType's that can be put in, it will automatically get the right dimension.
 * It is a wrapper around a pointer and redefines much used lapack and blas routines as memberfunctions
 */
class Vector
{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << vector << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << vector << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param vector_p de Vector you want to print
    */
   friend std::ostream &operator<<(std::ostream &output,const doci2DM::Vector &vector_p);

   public:

      //construct with as input a dimension
      Vector(int);

      //construct with as input a Matrix
      Vector(Matrix &);

      //copy constructor
      Vector(const Vector &);

      Vector(Vector &&);

      //destructor
      virtual ~Vector();

      //overload equality operator
      Vector &operator=(const Vector &);

      Vector &operator=(Vector &&);

      Vector &operator=(double);

      //overload += operator
      Vector &operator+=(const Vector &);

      //overload -= operator
      Vector &operator-=(const Vector &);

      Vector &daxpy(double alpha,const Vector &);

      Vector &operator/=(double);

      Vector &operator*=(double);

      void diagonalize(Matrix &);

      //easy to change the numbers
      double &operator[](int i);

      //easy to access the numbers
      double operator[](int i) const;

      void fill_Random();

      void fill_Random(int seed);

      //get the pointer to the vector
      double *gVector();

      //const version
      const double *gVector() const;

      int gn() const;

      double sum() const;

      double trace() const;

      double log_product() const;

      double ddot(const Vector &) const;

      void dscal(double alpha);

      double min() const;

      double max() const;

      void invert();

      void sqrt(int);

      void symmetrize();

      void L_map(const Vector &,const Vector &);

      Vector &mprod(const Vector &,const Vector &);

      void sort();

      void sep_pm(Vector &, Vector &);

   private:

      //!pointer of doubles, contains the numbers, the vector
      std::unique_ptr<double []> vector;

      //!dimension of the vector
      int n;
};

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
