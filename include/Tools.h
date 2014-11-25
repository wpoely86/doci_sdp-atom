#ifndef TOOLS_H
#define TOOLS_H

#include <string>

namespace doci2DM
{

class Tools
{
   public:
      static int getspDimension(std::string filename);

      static int getNumberOfParticles(std::string filename);

      static double getNuclearRepulEnergy(std::string filename);
};

}

#endif /* TOOLS_H */

/* vim: set ts=3 sw=3 expandtab :*/
