#ifndef POTENTIAL_REDUCTION_H
#define POTENTIAL_REDUCTION_H

#include "include.h"

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{

class PotentialReduction: public Method
{
   public:

      PotentialReduction(const CheMPS2::Hamiltonian &);

      PotentialReduction(const TPM &);

      PotentialReduction(const PotentialReduction &);

      PotentialReduction(PotentialReduction &&) = default;

      virtual ~PotentialReduction() = default;

      PotentialReduction* Clone() const;

      PotentialReduction* Move();

      PotentialReduction& operator=(PotentialReduction &&) = default;

      PotentialReduction& operator=(const PotentialReduction &);

      void BuildHam(const CheMPS2::Hamiltonian &);

      void BuildHam(const TPM &);

      unsigned int Run();

      double getFullEnergy() const;

      void set_target(double);

      void set_tolerance(double);

      void set_reduction(double);

      TPM& getRDM() const;

      TPM& getHam() const;

      Lineq& getLineq() const;

      double evalEnergy() const;

   private:

      std::unique_ptr<TPM> ham;

      std::unique_ptr<TPM> rdm;

      std::unique_ptr<Lineq> lineq;

      double nuclrep;

      double norm_ham;

      double target;

      double tolerance;

      double reductionfac;
};

}

#endif /* POTENTIAL_REDUCTION_H */

/* vim: set ts=3 sw=3 expandtab :*/
