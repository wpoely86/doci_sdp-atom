#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <memory>
#include <string>
#include <functional>

#include "Container.h"
#include "helpers.h"

namespace doci2DM
{

class SUP;
class Lineq;
class PHM;

class TPM: public Container
{
   friend std::ostream &operator<<(std::ostream &output,doci2DM::TPM &tpm);

   public:

      TPM(int L, int N);

      TPM(const TPM &) = default;

      TPM(TPM &&) = default;

      virtual ~TPM() = default;

      TPM& operator=(const TPM &) = default;

      TPM& operator=(TPM &&) = default;

      using Container::operator=;

      using Container::operator();

      double operator()(int a, int b, int c, int d) const;

      int gN() const;

      int gL() const;

      int gn() const;

      void HF_molecule(std::string filename);

      void ham(std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V);

      void WriteToFile(std::string filename) const;

      void ReadFromFile(std::string filename);

      double S_2() const;

      void unit();

      void init(const Lineq &);

      void Proj_Tr();

      void collaps(const SUP &, const Lineq &);

      void constr_grad(double t,const SUP &, const TPM &, const Lineq &);

      void Q(const TPM &);

      void Q(double a, double b, double c, const TPM &, bool=false);

      int solve(double t, const SUP &, TPM &, const Lineq &);

      void H(double t,const TPM &, const SUP &, const Lineq &);

      double line_search(double t, SUP &, const TPM &) const;

      double line_search(double t, const TPM &, const TPM &) const;

      void ReadFromFileFull(std::string filename);

      void Proj_E(const Lineq &, int option=0);

      std::vector<TPM> DOCI_constrains() const;

      std::vector<TPM> singlet_constrains() const;

      void S(const TPM &, bool=false);

      int InverseS(TPM &, const Lineq &);

      double getDiag(int, int) const;

      void G(const PHM &);

   private:

      void constr_lists(int L);

      //! number of particles
      int N;

      //! the size of the sp DOCI space (there are 2*L sp states)
      int L;

      //! dimension of the full TPM
      int n;

      //! table translating single particles indices to two particle indices
      static std::unique_ptr<helpers::matrix> s2t;

      //! table translating two particles indices to single particle indices
      static std::unique_ptr<helpers::matrix> t2s;
};

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
