#include "SUP.h"

using namespace doci2DM;

SUP::SUP(int L, int N)
{
   this->L = L;
   this->N = N;

   I.reset(new TPM(L,N));
#ifdef __Q_CON
   Q.reset(new TPM(L,N));
#endif

#ifdef __G_CON
   G.reset(new PHM(L,N));
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

#ifdef __G_CON
   G.reset(new PHM(*orig.G));
#endif
}

SUP::SUP(SUP &&orig)
{
   this->L = orig.L;
   this->N = orig.N;

   I = std::move(orig.I);
   Q = std::move(orig.Q);
   G = std::move(orig.G);
}

SUP& SUP::operator=(const SUP &orig)
{
   I.reset(new TPM(*orig.I));
#ifdef __Q_CON
   Q.reset(new TPM(*orig.Q));
#endif

#ifdef __G_CON
   G.reset(new PHM(*orig.G));
#endif

   return *this;
}

SUP& SUP::operator=(SUP &&orig)
{
   I = std::move(orig.I);
   Q = std::move(orig.Q);
   G = std::move(orig.G);

   return *this;
}

SUP& SUP::operator=(double a)
{
   (*I) = a;
#ifdef __Q_CON
   (*Q) = a;
#endif

#ifdef __G_CON
   (*G) = a;
#endif

   return *this;
}

SUP& SUP::operator+=(const SUP &orig)
{
   (*I) += (*orig.I);
#ifdef __Q_CON
   (*Q) += (*orig.Q);
#endif

#ifdef __G_CON
   (*G) += (*orig.G);
#endif

   return *this;
}

SUP& SUP::operator-=(const SUP &orig)
{
   (*I) -= (*orig.I);
#ifdef __Q_CON
   (*Q) -= (*orig.Q);
#endif

#ifdef __G_CON
   (*G) -= (*orig.G);
#endif

   return *this;
}

SUP& SUP::operator*=(double alpha)
{
   (*I) *= alpha;
#ifdef __Q_CON
   (*Q) *= alpha;
#endif

#ifdef __G_CON
   (*G) *= alpha;
#endif

   return *this;
}

SUP& SUP::operator/=(double alpha)
{
   (*I) /= alpha;
#ifdef __Q_CON
   (*Q) /= alpha;
#endif

#ifdef __G_CON
   (*G) /= alpha;
#endif

   return *this;
}

void SUP::dscal(double alpha)
{
   I->dscal(alpha);

#ifdef __Q_CON
   Q->dscal(alpha);
#endif

#ifdef __G_CON
   G->dscal(alpha);
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

PHM const & SUP::getG() const
{
   return *G;
}

PHM& SUP::getG()
{
   return *G;
}

void SUP::invert()
{
   I->invert();

#ifdef __Q_CON
   Q->invert();
#endif

#ifdef __G_CON
   G->invert();
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

#ifdef __G_CON
   G->G(tpm);
#endif
}

void SUP::sqrt(int option)
{
   I->sqrt(option);

#ifdef __Q_CON
   Q->sqrt(option);
#endif

#ifdef __G_CON
   G->sqrt(option);
#endif
}

void SUP::L_map(const SUP &A, const SUP &B)
{
   I->L_map(*A.I,*B.I);

#ifdef __Q_CON
   Q->L_map(*A.Q,*B.Q);
#endif

#ifdef __G_CON
   G->L_map(*A.G,*B.G);
#endif
}

int SUP::gnr() const
{
   int res = I->gnr();

#ifdef __Q_CON
   res += Q->gnr();
#endif

#ifdef __G_CON
   res += G->gnr();
#endif

   return res;
}

namespace doci2DM
{
   std::ostream &operator<<(std::ostream &output,doci2DM::SUP &sup)
   {
      output << "I block:" << std::endl;
      output << *sup.I << std::endl;

#ifdef __Q_CON
      output << "Q block:" << std::endl;
      output << *sup.Q << std::endl;
#endif

#ifdef __G_CON
      output << "G block:" << std::endl;
      output << *sup.G << std::endl;
#endif

      return output;
   }
}

double SUP::ddot(const SUP &x) const
{
   double result;

   result = I->ddot(*x.I);

#ifdef __Q_CON
   result += Q->ddot(*x.Q);
#endif

#ifdef __G_CON
   result += G->ddot(*x.G);
#endif

   return result;
}

void SUP::daxpy(double alpha, const SUP &y)
{
   I->daxpy(alpha, *y.I);
   
#ifdef __Q_CON
   Q->daxpy(alpha, *y.Q);
#endif

#ifdef __G_CON
   G->daxpy(alpha, *y.G);
#endif
}

void SUP::sep_pm(SUP &pos, SUP &neg)
{
   I->sep_pm(*pos.I, *neg.I);
   
#ifdef __Q_CON
   Q->sep_pm(*pos.Q, *neg.Q);
#endif

#ifdef __G_CON
   G->sep_pm(*pos.G, *neg.G);
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

/**
 * Write a SUP object to a HDF5 file
 * @param filename the name of the file to write to
 */
void SUP::WriteToFile(std::string filename) const
{
   hid_t       file_id, main_group_id, group_id;
   herr_t      status;

   file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   main_group_id = H5Gcreate(file_id, "SUP", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   group_id = H5Gcreate(main_group_id, "I", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   I->WriteToFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

#ifdef __Q_CON
   group_id = H5Gcreate(main_group_id, "Q", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   Q->WriteToFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);
#endif

#ifdef __G_CON
   group_id = H5Gcreate(main_group_id, "G", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   G->WriteToFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);
#endif

   status = H5Gclose(main_group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

/**
 * Read a SUP object from a HDF5 file. The SUP needs to be already created
 * with the correct dimensions. Will fill in the object on which it is called.
 * @param filename the file to read
 */
void SUP::ReadFromFile(std::string filename)
{
   hid_t       file_id, main_group_id, group_id;
   herr_t      status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   main_group_id = H5Gopen(file_id, "/SUP", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   group_id = H5Gopen(main_group_id, "I", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   I->ReadFromFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

#ifdef __Q_CON
   group_id = H5Gopen(main_group_id, "Q", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   Q->ReadFromFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);
#endif

#ifdef __G_CON
   group_id = H5Gopen(main_group_id, "G", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   G->ReadFromFile(group_id);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);
#endif

   status = H5Gclose(main_group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

/*  vim: set ts=3 sw=3 expandtab :*/
