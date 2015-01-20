/* Copyright (C) 2014  Ward Poelmans

   This file is part of Hubbard-GPU.

   Hubbard-GPU is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Hubbard-GPU is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Hubbard-GPU.  If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef HELPER_MATRIX_H
#define HELPER_MATRIX_H

#include <memory>
#include <complex>

namespace helpers {

/**
 * Helper class, wrapper around a double array. Has methods to get the number of rows
 * and columns. A simple matrix class.
 */
class matrix
{
    public:
        matrix();

        matrix(int n_, int m_);

        matrix(const matrix &orig);

        matrix(matrix &&orig);

        virtual ~matrix() { }

        matrix& operator=(const matrix &orig);

        matrix& operator=(double val);

        int getn() const;

        int getm() const;

        double operator()(int x,int y) const;

        double& operator()(int x,int y);

        double& operator[](int x);

        double operator[](int x) const;

        double* getpointer() const;

        matrix& prod(matrix const &A, matrix const &B);

        std::unique_ptr<double []> svd();

        std::unique_ptr<double []> sym_eig();

        void Print() const;

        double trace() const;

        void SaveToFile(std::string filename) const;

        void ReadFromFile(std::string filename) const;

    private:
        //!n by m array of double
        std::unique_ptr<double []> mat;
        //! number of rows
        int n;
        //! number of columns
        int m;
};

/**
 * Helper class, wrapper around a complex double array. Has methods to get the number of rows
 * and columns. A simple matrix class.
 */
class cmatrix
{
    public:
        cmatrix();

        cmatrix(int n_, int m_);

        cmatrix(const cmatrix &orig);

        cmatrix(cmatrix &&orig);

        virtual ~cmatrix() { }

        cmatrix& operator=(const cmatrix &orig);

        cmatrix& operator=(std::complex<double> val);

        int getn() const;

        int getm() const;

        std::complex<double> operator()(int x,int y) const;

        std::complex<double>& operator()(int x,int y);

        std::complex<double>& operator[](int x);

        std::complex<double> operator[](int x) const;

        std::complex<double>* getpointer() const;

        matrix& prod(cmatrix const &A, cmatrix const &B);

        void Print() const;

    private:
        //!n by m array of complex<double>
        std::unique_ptr<std::complex<double> []> mat;
        //! number of rows
        int n;
        //! number of columns
        int m;
};

template<typename T>
class tmatrix
{
    public:
        tmatrix();

        tmatrix(int n_, int m_);

        tmatrix(const tmatrix<T> &orig);

        tmatrix(tmatrix<T> &&orig);

        virtual ~tmatrix() { }

        tmatrix<T>& operator=(const tmatrix<T> &orig);

        tmatrix<T>& operator=(T val);

        unsigned int getn() const;

        unsigned int getm() const;

        T operator()(int x,int y) const;

        T& operator()(int x,int y);

        T& operator[](int x);

        T operator[](int x) const;

        T* getpointer() const;

        void Print() const;

    private:
        //!n by m array of double
        std::unique_ptr<T []> mat;
        //! number of rows
        unsigned int n;
        //! number of columns
        unsigned int m;
};

}

#endif /* HELPER_MATRIX_H */

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
