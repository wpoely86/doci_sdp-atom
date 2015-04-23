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

#include <assert.h>

#include "include.h"

using namespace doci2DM;

/**
 * standard constructor, only norm and DOCI constraints
 * @param L nr of levels
 * @param N nr of particles
 * @param partial_trace use constraints for the block and vector seperatly
 */
Lineq::Lineq(int L,int N, bool partial_trace)
{
   this->L = L;
   this->N = N;

   if(partial_trace)
   {
      TPM trace_constr(L,N);
      trace_constr = 0;

      for(int i=0;i<trace_constr.gdimMatrix(0);i++)
         trace_constr(0,i,i) = 1;

      E.push_back(trace_constr);
      e.push_back(N/2.0);

      trace_constr.getMatrix(0) = 0;
      for(int i=0;i<trace_constr.gdimVector(0);i++)
         trace_constr(0,i) = 1;

      E.push_back(std::move(trace_constr));
      e.push_back(N*(N/2.0-1)); // keep the fourfold degenaracy in mind
   } else
   {
      TPM trace_constr(L,N);
      trace_constr.unit();

      // make room for 1 trace constrains and L DOCI constrains
      e.resize(1+L, 0);

      // first trace on Block
      E.push_back(std::move(trace_constr));
      e[0] = N*(N-1)/2.0;

      auto doci_constr = trace_constr.DOCI_constrains();
      E.insert(E.end(), std::make_move_iterator(doci_constr.begin()), std::make_move_iterator(doci_constr.end()));

      assert(E.size() == 1+L);
   }

   orthogonalize();//speaks for itself, doesn't it?

   //construct the u^0 SUP matrices
   constr_u_0();

   orthogonalize_u_0();
}

/**
 * @return nr of particles
 */
int Lineq::gN() const
{
   return N;
}

/**
 * @return nr of constraints
 */
int Lineq::gnr() const
{
   return E.size();
}

/**
 * @return nr of sp orbitals
 */
int Lineq::gL() const
{
   return L;
}

/**
 * access to the individual constraint TPM's
 * @param i the index
 * @return The E TPM on index i.
 */
const TPM &Lineq::gE(int i) const
{
   return E[i];
}

/**
 * access to the individual constraint values
 * @param i the index
 * @return the e values on index i: e[i] or something
 */
double Lineq::ge(int i) const
{
   return e[i];
}

/**
 * access to the individual orthogonalized constraint TPM's: private function
 * @param i the index
 * @return The E_ortho TPM on index i.
 */
const TPM &Lineq::gE_ortho(int i) const
{
   return E_ortho[i];
}

/**
 * access to the individual orthogonalized constraint values: private function
 * @param i the index
 * @return the e values on index i: e_ortho[i] or something
 */
double Lineq::ge_ortho(int i) const
{
   return e_ortho[i];
}

namespace doci2DM
{
   std::ostream &operator<<(std::ostream &output,doci2DM::Lineq &lineq_p)
   {
      output << "first print the constraint matrices:";
      output << std::endl;
      output << std::endl;

      for(int i = 0;i < lineq_p.gnr();++i)
      {
         output << "constraint nr :" << i << std::endl;
         output << std::endl;

         output << lineq_p.gE(i);
      }

      output << std::endl;
      output << std::endl;
      output << "the desired values are:" << std::endl;
      output << std::endl;

      for(int i = 0;i < lineq_p.gnr();++i)
         output << i << "\t" << lineq_p.ge(i) << std::endl;

      return output;
   }
}

/**
 * orthogonalize the constraints, will take E and e, and construct E_ortho and e_ortho with them.
 */
void Lineq::orthogonalize()
{
   //construct the overlapmatrix of the E's
   Matrix S(E.size());

   for(int i = 0;i < gnr();++i)
      for(int j = i;j < gnr();++j)
         S(i,j) = S(j,i) = E[i].ddot(E[j]);

   //take the inverse square root
   S.sqrt(-1);

   // make room
   E_ortho.reserve(gnr());
   e_ortho.reserve(gnr());
   // throw away old stuff (leave capacity unchanged)
   E_ortho.clear();
   e_ortho.clear();

   //make the orthogonal ones:
   for(int i=0;i<gnr();++i)
   {
      TPM ortho_constr(L,N);
      ortho_constr = 0;

      double tmp_e_ortho = 0;

      for(int j = 0;j<gnr();++j)
      {
         ortho_constr.daxpy(S(i,j), E[j]);
         tmp_e_ortho += S(i,j) * e[j];
      }

      E_ortho.push_back(std::move(ortho_constr));
      e_ortho.push_back(tmp_e_ortho);
   }

   assert(E_ortho.size() == E.size());
}

/**
 * access to the u_0 matrices
 * @param i the index of the specific u_0 matrix you are interested in.
 * @return u_0[i]
 */
const SUP &Lineq::gu_0(int i) const
{
   return u_0[i];
}

/**
 * access to the orthogonalized u_0 matrices
 * @param i the index of the specific u_0_ortho matrix you are interested in.
 * @return u_0_ortho[i]
 */
const SUP &Lineq::gu_0_ortho(int i) const
{
   return u_0_ortho[i];
}

/**
 * construct's the u_0 matrices with the input provided in the constructor. Again reusable code for the different constructers.
 * and again private because I say so.
 */
void Lineq::constr_u_0()
{
   u_0.reserve(gnr());

   for(int i = 0;i < gnr();++i)
   {
      SUP tmp(L,N);

      tmp.fill(E_ortho[i]);

      u_0.push_back(std::move(tmp));
   }

   assert(u_0.size() == E.size());
}

/**
 * orthogonalize the u_0 matrices, will take the u_0's and construct the u_0_ortho's.
 */
void Lineq::orthogonalize_u_0()
{
   //construct the overlapmatrix of the E's
   Matrix S(u_0.size());

   for(int i = 0;i < gnr();++i)
      for(int j = i;j < gnr();++j)
         S(i,j) = S(j,i) = u_0[i].ddot(u_0[j]);

   //take the inverse square root
   S.sqrt(-1);

   // make room
   u_0_ortho.reserve(gnr());
   // throw away old stuff (leave capacity unchanged)
   u_0_ortho.clear();

   //make the orthogonal ones:
   for(int i=0;i<gnr();++i)
   {
      SUP ortho_constr(L,N);
      ortho_constr = 0;

      for(int j = 0;j<gnr();++j)
         ortho_constr.daxpy(S(i,j), u_0[j]);

      u_0_ortho.push_back(std::move(ortho_constr));
   }

   assert(u_0_ortho.size() == E.size());
//   //first allocate the u_0_ortho matrices:
//   u_0_ortho = new SUP * [nr];
//
//   for(int i = 0;i < nr;++i)
//      u_0_ortho[i] = new SUP(M,N);
//
//   //construct the overlapmatrix of the u_0's
//   Matrix S(nr);
//
//   for(int i = 0;i < nr;++i)
//      for(int j = i;j < nr;++j)
//         S(i,j) = u_0[i]->ddot(*u_0[j]);
//
//   S.symmetrize();
//
//   //take the inverse square root
//   S.sqrt(-1);
//
//   //make the orthogonal ones:
//   for(int i = 0;i < nr;++i){
//
//      *u_0_ortho[i] = 0;//init
//
//      for(int j = 0;j < nr;++j)
//         u_0_ortho[i]->daxpy(S(i,j),*u_0[j]);
//
//   }
}

void Lineq::check(const TPM &tpm) const
{
   for(int i = 0;i < gnr();++i)
   {
       std::cout << "constrain " << i << " gives " << tpm.ddot(gE(i)) << " = " << ge(i);
       double tmp = fabs(tpm.ddot(gE(i)) - ge(i));
       std::cout << ((tmp>1e-10) ? "\tFAILED" : "") << std::endl;
   }
}

/*  vim: set ts=3 sw=3 expandtab :*/
