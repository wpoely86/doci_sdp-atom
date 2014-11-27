#ifndef METHOD_H
#define METHOD_H

#include <memory>

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{
class TPM;

class Method
{
   public:

      virtual ~Method() = default;

      virtual void Run() = 0;

      virtual void BuildHam(const CheMPS2::Hamiltonian &) = 0;

      virtual Method* Clone() const = 0;

      virtual Method* Move() = 0;

      virtual TPM& getRDM() const = 0;

      double getEnergy() const { return energy; }

   protected:

      int L;

      int N;

      double energy;
};

}

#endif /* METHOD_H */

/* vim: set ts=3 sw=3 expandtab :*/
