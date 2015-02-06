#ifndef TOOLS_H
#define TOOLS_H

#include <string>

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{
   class TPM;

class Tools
{
   public:
      static int getspDimension(std::string filename);

      static int getNumberOfParticles(std::string filename);

      static double getNuclearRepulEnergy(std::string filename);

      static void scan_all(const TPM &rdm, const CheMPS2::Hamiltonian &ham);
};

}

#endif /* TOOLS_H */

/* vim: set ts=3 sw=3 expandtab :*/
