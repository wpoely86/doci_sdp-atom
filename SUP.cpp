#include "SUP.h"


SUP::SUP(int L, int N)
{
   this->L = L;
   this->N = N;

   I.reset(new TPM(L,N));
#ifdef __Q_CON
   Q.reset(new TPM(L,N));
#endif
}

SUP::SUP(const SUP &orig)
{
   this->L = orig.L;
   this->N = orig.N;

   I.reset(new TPM(*orig.I));
#ifdef __Q_CON
   Q.reset(new TPM(*orig.Q));
#endif
}

SUP::SUP(SUP &&orig)
{
   this->L = orig.L;
   this->N = orig.N;

   I = std::move(orig.I);
   Q = std::move(orig.Q);
}

SUP& SUP::operator=(const SUP &orig)
{
   I.reset(new TPM(*orig.I));
#ifdef __Q_CON
   Q.reset(new TPM(*orig.Q));
#endif

   return *this;
}

SUP& SUP::operator=(SUP &&orig)
{
   I = std::move(orig.I);
   Q = std::move(orig.Q);

   return *this;
}

SUP& SUP::operator=(double a)
{
   (*I) = a;
#ifdef __Q_CON
   (*Q) = a;
#endif

   return *this;
}

SUP& SUP::operator+=(const SUP &orig)
{
   (*I) += (*orig.I);
#ifdef __Q_CON
   (*Q) += (*orig.Q);
#endif

   return *this;
}

SUP& SUP::operator-=(const SUP &orig)
{
   (*I) -= (*orig.I);
#ifdef __Q_CON
   (*Q) -= (*orig.Q);
#endif

   return *this;
}

SUP& SUP::operator*=(double alpha)
{
   (*I) *= alpha;
#ifdef __Q_CON
   (*Q) *= alpha;
#endif

   return *this;
}

SUP& SUP::operator/=(double alpha)
{
   (*I) /= alpha;
#ifdef __Q_CON
   (*Q) /= alpha;
#endif

   return *this;
}

void SUP::dscal(double alpha)
{
   I->dscal(alpha);

#ifdef __Q_CON
   Q->dscal(alpha);
#endif
}

int SUP::gN() const
{
   return N;
}

int SUP::gL() const
{
   return L;
}

TPM const & SUP::getI() const
{
   return *I;
}

TPM& SUP::getI()
{
   return *I;
}

TPM const & SUP::getQ() const
{
   return *Q;
}

TPM& SUP::getQ()
{
   return *Q;
}

void SUP::invert()
{
   I->invert();

#ifdef __Q_CON
   Q->invert();
#endif
}

/**
 * Fill the SUP object with the tpm object
 * @param tpm the TPM to use
 */
void SUP::fill(const TPM &tpm)
{
   (*I) = tpm;

#ifdef __Q_CON
   Q->Q(tpm);
#endif
}

void SUP::sqrt(int option)
{
   I->sqrt(option);

#ifdef __Q_CON
   Q->sqrt(option);
#endif
}

void SUP::L_map(const SUP &A, const SUP &B)
{
   I->L_map(*A.I,*B.I);

#ifdef __Q_CON
   Q->L_map(*A.Q,*B.Q);
#endif
}

int SUP::gnr() const
{
   int res = I->gnr();

#ifdef __Q_CON
   res += Q->gnr();
#endif

   return res;
}

std::ostream &operator<<(std::ostream &output,SUP &sup)
{
   output << "I block:" << std::endl;
   output << *sup.I << std::endl;

#ifdef __Q_CON
   output << "Q block:" << std::endl;
   output << *sup.Q << std::endl;
#endif

   return output;
}

double SUP::ddot(const SUP &x) const
{
   double result;

   result = I->ddot(*x.I);

#ifdef __Q_CON
   result += Q->ddot(*x.Q);
#endif

   return result;
}

void SUP::daxpy(double alpha, const SUP &y)
{
   I->daxpy(alpha, *y.I);
   
#ifdef __Q_CON
   Q->daxpy(alpha, *y.Q);
#endif
}

void SUP::sep_pm(SUP &pos, SUP &neg)
{
   I->sep_pm(*pos.I, *neg.I);
   
#ifdef __Q_CON
   Q->sep_pm(*pos.Q, *neg.Q);
#endif
}

/**
 * Initialization of the SUP matrix S, is just u^0: see primal_dual.pdf for more information
 */
void SUP::init_S(const Lineq &lineq)
{
   *this = lineq.gu_0(0);

   this->dscal(lineq.ge_ortho(0));

   for(int i = 1;i < lineq.gnr();++i)
      this->daxpy(lineq.ge_ortho(i),lineq.gu_0(i));
}

/*  vim: set ts=3 sw=3 expandtab :*/
