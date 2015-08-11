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

#ifndef LINEQ_H
#define LINEQ_H

#include <iostream>
#include <vector>

#include "TPM.h"
#include "SUP.h"

namespace doci2DM
{

/**
 * @author Brecht Verstichel
 * @date 09-09-2010\n\n
 * This is a class written for handling linear equality constraints in the program. Contains a collection of the matrices
 * on which the tpm has to have a certain projection, these values are stored in a pointer of doubles.
 */

class Lineq
{
   /**
    * Output stream operator overloaded: prints all the constraint matrices with a whitespace between them.
    * @param output The stream to which you are writing (e.g. cout)
    * @param lineq_p The Lineq you want to print
    */
   friend std::ostream &operator<<(std::ostream &output,Lineq &lineq_p);

   public:

      Lineq(int L,int N, bool=false);

      virtual ~Lineq() = default;

      int gN() const;

      int gnr() const;

      int gL() const;

      const TPM &gE(int) const;

      double ge(int) const;

      const TPM &gE_ortho(int) const;

      double ge_ortho(int) const;

      const SUP &gu_0(int) const;

      const SUP &gu_0_ortho(int) const;

      void check(const TPM &tpm) const;

      void orthogonalize();

   private:

      void constr_u_0();

      void orthogonalize_u_0();

      //!double pointer to TPM object, will contain the linear equality constraints
      std::vector<TPM> E;
      //!pointer of doubles, will contain the values of the projections. (the desired equalities)
      std::vector<double> e;
      
      //!orthogonalized constraints, these will be hidden from the public.
      std::vector<TPM> E_ortho;
      //!the values accompanying the orthogonalized constraints
      std::vector<double> e_ortho;

      //!will contain the matrices that span u^0 space
      std::vector<SUP> u_0;

      //!will contain the orthogonalized matrices that span u^0 space
      std::vector<SUP> u_0_ortho;

      //!nr of particles
      int N;

      //!nr of sp orbs
      int L;
};

}

#endif /* LINEQ_H */

/* vim: set ts=3 sw=3 expandtab :*/
