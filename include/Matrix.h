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

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>
#include <memory>

namespace doci2DM
{

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written for symmetric matrices. It is a wrapper around a double pointer and
 * redefines much used lapack and blas routines as memberfunctions
 */

class Vector;

class Matrix
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de Matrix you want to print
    */
   friend std::ostream &operator<<(std::ostream &,const doci2DM::Matrix &);

   public:

      //constructor
      Matrix(int n);

      //copy constructor
      Matrix(const Matrix &);

      Matrix(Matrix &&);

      //destructor
      virtual ~Matrix();

      //overload equality operator
      Matrix &operator=(const Matrix &);

      Matrix &operator=(Matrix &&);

      Matrix &operator=(double );

      //overload += operator
      Matrix &operator+=(const Matrix &);

      //overload -= operator
      Matrix &operator-=(const Matrix &);

      Matrix &daxpy(double alpha,const Matrix &);

      Matrix &operator*=(double);

      Matrix &operator/=(double);

      Matrix &mprod(const Matrix &,const Matrix &);

      //easy to change the numbers
      double &operator()(int i,int j);

      //easy to access the numbers
      double operator()(int i,int j) const;

      //get the pointer to the matrix
      double *gMatrix();

      const double *gMatrix() const;

      int gn() const;

      double trace() const;

      Vector diagonalize();

      Vector diagonalize_2x2();

      double ddot(const Matrix &) const;

      void invert();

      void invert_2x2();

      void dscal(double alpha);

      void fill_Random();

      void fill_Random(int seed);

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void sqrt_2x2(int option);

      void mdiag(const Vector &);

      void L_map(const Matrix &,const Matrix &);

      void L_map_2x2(const Matrix &,const Matrix &);

      void symmetrize();

      void SaveRawToFile(const std::string) const;

      void sep_pm(Matrix &,Matrix &);

      void sep_pm_2x2(Matrix &,Matrix &);

      void unit();

   private:

      //!pointer of doubles, contains the numbers, the matrix
      std::unique_ptr<double []> matrix;

      //!dimension of the matrix
      int n;
};

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
