#ifndef EIG_H
#define EIG_H

#include <iostream>

#include "include.h"

namespace doci2DM
{

class EIG: public BlockVector
{
   public:

      EIG(BlockMatrix &);

      EIG(Container &);

      EIG(SUP &);

      virtual ~EIG() = default;

      using BlockVector::operator=;

      using BlockVector::operator();

      using BlockVector::operator[];

      double min() const;

      double max() const;

      double lsfunc(double) const;
};

}

#endif /* EIG_H */

/* vim: set ts=3 sw=3 expandtab :*/
