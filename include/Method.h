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

      Method() { do_output = true; }

      virtual ~Method() = default;

      virtual unsigned int Run() = 0;

      virtual void BuildHam(const CheMPS2::Hamiltonian &) = 0;

      virtual Method* Clone() const = 0;

      virtual Method* Move() = 0;

      virtual TPM& getRDM() const = 0;

      double getEnergy() const { return energy; }

      virtual void set_output(bool out) { do_output = out; } 

      virtual void set_outfile(std::string filename) { outfile = filename; }

   protected:

      int L;

      int N;

      double energy;

      bool do_output;

      std::string outfile;
};

}

#endif /* METHOD_H */

/* vim: set ts=3 sw=3 expandtab :*/
