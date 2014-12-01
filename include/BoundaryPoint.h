#ifndef BOUNDARY_POINT_H
#define BOUNDARY_POINT_H

#include "include.h"

namespace CheMPS2 { class Hamiltonian; }

namespace doci2DM
{

class BoundaryPoint: public Method
{
   public:

      BoundaryPoint(const CheMPS2::Hamiltonian &);

      BoundaryPoint(const BoundaryPoint &);

      BoundaryPoint(BoundaryPoint &&) = default;

      virtual ~BoundaryPoint() = default;

      BoundaryPoint& operator=(const BoundaryPoint &);

      BoundaryPoint& operator=(BoundaryPoint &&) = default;

      BoundaryPoint* Clone() const;

      BoundaryPoint* Move();

      void BuildHam(const CheMPS2::Hamiltonian &);

      unsigned int Run();

      double getFullEnergy() const;

      void set_tol_PD(double);

      void set_tol_en(double);

      void set_mazzy(double);

      void set_sigma(double);

      void set_max_iter(int);

      SUP& getX() const;

      SUP& getZ() const;

      Lineq& getLineq() const;

      TPM& getRDM() const;

   private:

      std::unique_ptr<TPM> ham;

      std::unique_ptr<SUP> X;

      std::unique_ptr<SUP> Z;

      std::unique_ptr<Lineq> lineq;

      double nuclrep;

      double tol_PD, tol_en;

      double mazzy;

      double sigma;

      int max_iter;

      unsigned int avg_iters;

      unsigned int iters;
      unsigned int runs;
};

}

#endif /* BOUNDARY_POINT_H */

/* vim: set ts=3 sw=3 expandtab :*/
