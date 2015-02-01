#include <cstdio>
#include <sstream>
#include <assert.h>
#include <hdf5.h>

#include "include.h"

using namespace doci2DM;

// default empty
std::unique_ptr<helpers::tmatrix<unsigned int>> TPM::s2t = nullptr;
std::unique_ptr<helpers::tmatrix<unsigned int>> TPM::t2s = nullptr;

/**
 * Create a object in the Two Particle space
 * @param L the number of levels (the sp space has size of 2*L)
 * @param N the number of particles
 */
TPM::TPM(int L, int N): Container(1,1)
{
   this->N = N;
   this->L = L;
   n = L*(2*L-1);

   // the LxL block with degen = 1
   setMatrixDim(0, L, 1);
   // the L/2*(L-1) vector with degen = 4
   setVectorDim(0, (L*(L-1))/2, 4);

   if(!s2t || !t2s)
      constr_lists(L);
}

void TPM::constr_lists(int L)
{
   int M = 2*L;
   int n_tp = M*(M-1)/2;

   s2t.reset(new helpers::tmatrix<unsigned int>(M,M));
   (*s2t) = -1; // if you use something you shouldn't, this will case havoc

   t2s.reset(new helpers::tmatrix<unsigned int>(n_tp,2));
   (*t2s) = -1; // if you use something you shouldn't, this will case havoc

   int tel = 0;

   // a \bar a
   for(int a=0;a<L;a++)
      (*s2t)(a,a+L) = (*s2t)(a+L,a) = tel++;

   // a b
   for(int a=0;a<L;a++)
      for(int b=a+1;b<L;b++)
         (*s2t)(a,b) = (*s2t)(b,a) = tel++;

   // \bar a \bar b
   for(int a=L;a<M;a++)
      for(int b=a+1;b<M;b++)
         (*s2t)(a,b) = (*s2t)(b,a) = tel++;

   // a \bar b ; a \bar b
   for(int a=0;a<L;a++)
      for(int b=L+a+1;b<M;b++)
         if(a%L!=b%L)
            (*s2t)(a,b) = (*s2t)(b,a) = tel++;

   // \bar a b ; \bar a b
   for(int a=L;a<M;a++)
      for(int b=a%L+1;b<L;b++)
         if(a%L!=b%L)
            (*s2t)(a,b) = (*s2t)(b,a) = tel++;

   assert(tel == n_tp);

   for(int a=0;a<M;a++)
      for(int b=a+1;b<M;b++)
      {
         (*t2s)((*s2t)(a,b),0) = a;
         (*t2s)((*s2t)(a,b),1) = b;
      }
}

int TPM::gN() const
{
   return N;
}

int TPM::gL() const
{
   return L;
}

int TPM::gn() const
{
   return n;
}

/**
 * General access operator to the rdm. Returns 0 if the sp indices
 * are not in DOCI space
 * @param a first sp index
 * @param b second sp index
 * @param c thirth sp index
 * @param d fourth sp index
 * @returns the value of the rdm with sp index a,b|c,d
 */
double TPM::operator()(int a, int b, int c, int d) const
{
   if(a==b || c==d)
      return 0;

   int sign = 1;

   if(a > b)
      sign *= -1;

   if(c > d)
      sign *= -1;

   int i = (*s2t)(a,b);
   int j = (*s2t)(c,d);

   if(i<L && j<L)
      return sign * (*this)(0,i,j);
   else if(i==j)
      return sign * (*this)(0,(i-L)%gdimVector(0));
   else
      return 0;
}

namespace doci2DM {
   std::ostream &operator<<(std::ostream &output,doci2DM::TPM &tpm)
   {
      output << "Block: " << std::endl;
      for(int i=0;i<tpm.L;i++)
         for(int j=i;j<tpm.L;j++)
            output << i << "\t" << j << "\t|\t" << (*tpm.t2s)(i,0) << "  " <<  (*tpm.t2s)(i,1) << " ; " <<  (*tpm.t2s)(j,0) << "  " <<  (*tpm.t2s)(j,1) << "\t\t" << tpm(0,i,j) << std::endl;

      output << std::endl;

      output << "Vector (4x): " << std::endl;
      for(int i=0;i<tpm.getVector(0).gn();i++)
         output << i << "\t|\t" << (*tpm.t2s)(tpm.L+i,0) << "  " << (*tpm.t2s)(tpm.L+i,1) << "\t\t" << tpm(0,i) << std::endl;

      return output;
   }
}

/**
 * Build the reduced hamiltonian using 2 function objects: one for the one particles integrals,
 * one for the two particle integrals. The functions should accept the spatial orbitals
 * as argument.
 * @param T a function that accepts 2 orbitals a and b and returns the <a|T|b> matrix element
 * @param V a function that accepts 4 orbitals a, b, c and d and returns the <ab|V|cd> matrix element
 */
void TPM::ham(std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V)
{
   // make our life easier
   auto calc_elem = [this,&T,&V] (int i, int j) {
      int a = (*t2s)(i,0);
      int b = (*t2s)(i,1);
      int c = (*t2s)(j,0);
      int d = (*t2s)(j,1);

      int a_ = a % L;
      int b_ = b % L;
      int c_ = c % L;
      int d_ = d % L;

      double result = 0;

      // sp terms
      if(i==j)
         result += (T(a_,a_) + T(b_,b_))/(N - 1.0);

      // tp terms:

      // a \bar a ; b \bar b
      if(b==(a+L) && d==(c+L))
         result += V(a_,b_,c_,d_);

      // a b ; a b
      if(i==j && a<L && b<L && a!=b)
         result += V(a_,b_,c_,d_) - V(a_,b_,d_,c_);

      // \bar a \bar b ; \bar a \bar b
      if(i==j && a>=L && b>=L && a!=b)
         result += V(a_,b_,c_,d_) - V(a_,b_,d_,c_);

      // a \bar b ; a \bar b
      if(i==j && a<L && b>=L && a%L!=b%L)
         result += V(a_,b_,c_,d_);

      // \bar a b ; \bar a b
      if(i==j && a>=L && b<L && a%L!=b%L)
         result += V(a_,b_,c_,d_);

      return result;
   };


   auto& hamB = getMatrix(0);

   for(int i=0;i<L;++i)
      for(int j=i;j<L;++j)
         hamB(i,j) = hamB(j,i) = calc_elem(i,j);

   auto& hamV = getVector(0);

   for(int i=0;i<hamV.gn();i++)
      // keep in mind that the degen of the vector is 4. We need prefactor of 2, so
      // we end up with 0.5
      hamV[i] = 0.5*calc_elem(L+i,L+i) + 0.5*calc_elem(L*L+i,L*L+i);
}

void TPM::HF_molecule(std::string filename)
{
   hid_t       file_id, group_id, dataset_id;
   herr_t      status;

   // open file
   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/integrals", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   dataset_id = H5Dopen(group_id, "OEI", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   Matrix OEI(L);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, OEI.gMatrix());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   Matrix TEI(L*L);

   dataset_id = H5Dopen(group_id, "TEI", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, TEI.gMatrix());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

//   printf("0\t0\t0\t0\t0\n");
//   for(int a=0;a<L;a++)
//      for(int b=0;b<L;b++)
//         printf("%20.15f\t%d\t%d\t0\t0\n", OEI(a,b), a+1,b+1);
//
//   for(int a=0;a<L;a++)
//      for(int b=0;b<L;b++)
//         for(int c=0;c<L;c++)
//            for(int d=0;d<L;d++)
//               printf("%20.15f\t%d\t%d\t%d\t%d\n", TEI(a*L+b,c*L+d), a+1,c+1,b+1,d+1);

   std::function<double(int,int)> getT = [&OEI] (int a, int b) -> double { return OEI(a,b); };
   std::function<double(int,int,int,int)> getV = [&TEI,this] (int a, int b, int c, int d) -> double { return TEI(a*L + b,c*L + d); };

   ham(getT,getV);
}


void TPM::WriteToFile(std::string filename) const
{
   hid_t       file_id, group_id, dataset_id, attribute_id, dataspace_id;
   hsize_t     dims;
   herr_t      status;

   file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   group_id = H5Gcreate(file_id, "/RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);



   hsize_t dimblock = getMatrix(0).gn() * getMatrix(0).gn();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Block", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   double *data = const_cast<Matrix &>(getMatrix(0)).gMatrix();

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dimblock = getVector(0).gn();

   dataspace_id = H5Screate_simple(1, &dimblock, NULL);

   dataset_id = H5Dcreate(group_id, "Vector", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   data = const_cast<Vector &>(getVector(0)).gVector();

   status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dataspace_id = H5Screate(H5S_SCALAR);

   attribute_id = H5Acreate (group_id, "L", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &L );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   attribute_id = H5Acreate (group_id, "N", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &N );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);


   dims = 6;
   dataspace_id = H5Screate_simple(1, &dims, NULL);

   int typeofcalculation[6];

   typeofcalculation[0] = 1; // P

#ifdef __Q_CON
   typeofcalculation[1] = 1; // Q
#else
   typeofcalculation[1] = 0; // Q
#endif

#ifdef __G_CON
   typeofcalculation[2] = 1; // G
#else
   typeofcalculation[2] = 0; // Q
#endif

#ifdef __T1_CON
   typeofcalculation[3] = 1; // T1
#else
   typeofcalculation[3] = 0; // Q
#endif

#ifdef __T2_CON
   typeofcalculation[4] = 1; // T2
#else
   typeofcalculation[4] = 0; // Q
#endif

#ifdef __T2P_CON
   typeofcalculation[5] = 1; // T2P
#else
   typeofcalculation[5] = 0; // Q
#endif

   attribute_id = H5Acreate (group_id, "Type", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite (attribute_id, H5T_NATIVE_INT, &typeofcalculation[0] );
   HDF5_STATUS_CHECK(status);

   status = H5Aclose(attribute_id);
   HDF5_STATUS_CHECK(status);

   status = H5Sclose(dataspace_id);
   HDF5_STATUS_CHECK(status);

   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

void TPM::ReadFromFile(std::string filename)
{
   hid_t       file_id, group_id, dataset_id;
   herr_t      status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   group_id = H5Gopen(file_id, "/RDM", H5P_DEFAULT);
   HDF5_STATUS_CHECK(group_id);

   dataset_id = H5Dopen(group_id, "Block", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, getMatrix(0).gMatrix());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   dataset_id = H5Dopen(group_id, "Vector", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, getVector(0).gVector());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);


   status = H5Gclose(group_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);
}

/**
 * @todo Optimize to DOCI form
 * @returns S^2 expectation value
 */
double TPM::S_2() const
{
   double spin = 0;

   for(int i = 0;i<n;++i)
   {
      int a = (*t2s)(i,0);
      int b = (*t2s)(i,1);

      const double s_a = ( 1.0 - 2 * (a / L) )/2.0;
      const double s_b = ( 1.0 - 2 * (b / L) )/2.0;

      spin += ( (1 + s_a*s_a + s_b*s_b)/(N - 1.0) + 2*s_a*s_b ) * (*this)(a,b,a,b);
   }

   //then the off diagonal elements: a and b are sp indices
   for(int a = 0;a < L;++a)
      for(int b = 0;b < L;++b)
         spin += (*this)(a,L+b,a+L,b);

   return spin;
}

/**
 * Set this object to the unit matrix
 */
void TPM::unit()
{
   (*this) = 0;

   for(int i=0;i<gdimMatrix(0);i++)
      (*this)(0,i,i) = 1;

   for(int i=0;i<gdimVector(0);i++)
      (*this)(0,i) = 1;
}

/**
 * Construct a TPM that is OK with the constrains
 * @param lineq the constrains used in the optimalization
 */
void TPM::init(const Lineq &lineq)
{
   (*this) = 0;

   for(int i=0;i<lineq.gnr();i++)
      daxpy(lineq.ge_ortho(i), lineq.gE_ortho(i));
}

void TPM::Proj_Tr()
{
   assert(0 && "You should not be here!");

   double trace = getMatrix(0).trace()/L;

   for(int i=0;i<L;i++)
      (*this)(0,i,i) -= trace;

   trace = getVector(0).trace()*4.0/(n-L);

   for(int i=0;i<gdimVector(0);i++)
      (*this)(0,i) -= trace;
}

/**
 * Collaps the full SUP to the TPM space:
 * *this = S.I + Q(S.Q) + G^Down(S.G) + ...
 * stores the result in *this
 * @param S The SUP to collaps
 */
void TPM::collaps(const SUP &S, const Lineq &lineq)
{
   (*this) = S.getI();

#ifdef __Q_CON
   TPM hulp(L,N);
   hulp.Q(S.getQ());

   (*this) += hulp;
#endif

#ifdef __G_CON
   hulp.G(S.getG());

   *this += hulp;
#endif

   Proj_E(lineq);
}

/**
 * Construct the -gradient, store in *this
 * @param t the barrier height
 * @param S the SUP to use
 * @param ham the hamiltonian
 */
void TPM::constr_grad(double t,const SUP &S, const TPM &ham, const Lineq &lineq)
{
   collaps(S, lineq);

   (*this) *= t;

   (*this) -= ham;

   Proj_E(lineq);
}

/**
 * The Q image: stores the result in *this
 * @param tpm The input matrix
 */
void TPM::Q(const TPM &tpm)
{
   (*this) = tpm;

   double tmp = trace() / (N*(N-1)/2.0);

   // the matrix inequality
   for(int i=0;i<L;i++)
      (*this)(0,i,i) += tmp - 2 * tpm(0,i,i);

   SPM spm(L,N);
   spm.bar2(1.0/(N/2.0-1.0),tpm);

   // the linear inequality
   for(int i=0;i<getVector(0).gn();i++)
   {
      int a = (*t2s)(L+i,0);
      int b = (*t2s)(L+i,1);

      (*this)(0,i) += tmp - spm(0,a) - spm(0,b);
   }

//   Q(1, 2.0/(N*(N-1)), 1.0/(N-1.0),tpm);
}

void TPM::Q(double alpha, double beta, double gamma, const TPM &tpm, bool inverse)
{
   if(inverse)
   {
      beta = (beta*alpha + beta*gamma*2*L - 2.0*gamma*gamma)/(alpha * (gamma*(2*L - 2.0) -  alpha) * ( alpha + beta*2*L*(2*L - 1.0) - 2.0*gamma*(2*L - 1.0) ) );
      gamma = gamma/(alpha*(gamma*((2*L) - 2.0) - alpha));
      alpha = 1.0/alpha;
   }

   SPM spm(L,N);
   spm.bar3(gamma, tpm);

   (*this) = tpm;
   (*this) *= alpha;

   double tmp = beta * trace();

   // the matrix inequality
   for(int i=0;i<L;i++)
      (*this)(0,i,i) += tmp - 2 * spm(0,i);

   // the linear inequality
   for(int i=0;i<getVector(0).gn();i++)
   {
      int a = (*t2s)(L+i,0);
      int b = (*t2s)(L+i,1);

      (*this)(0,i) += tmp - spm(0,a) - spm(0,b);
   }
}

int TPM::solve(double t, const SUP &S, TPM &grad, const Lineq &lineq, int max_iters)
{
   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   TPM r(grad);

   //norm van het residu
   double rr = r.ddot(r);

   double rr_old, ward;

   TPM Hb(L,N);

   while(rr > 1.0e-9)
   { 
      Hb.H(t,grad,S, lineq);

      ward = rr/grad.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward, grad);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      grad.dscal(rr/rr_old);

      grad += r;

      ++iter;

      if(iter>max_iters)
         break;
   }

   return iter;
}

/**
 * Calculate the hessian
 * @param t barrier height
 * @param delta the current delta in the TPM space
 * @param S the full filled SUP object
 */
void TPM::H(double t,const TPM &delta,const SUP &S, const Lineq &lineq)
{
   L_map(S.getI(),delta);

#ifdef __Q_CON

   TPM tmp1(L,N);
   tmp1.Q(delta);

   TPM tmp2(L,N);
   tmp2.L_map(S.getQ(), tmp1);

   tmp1.Q(tmp2);

   (*this) += tmp1;

#endif

#ifdef __G_CON
   //hulpje voor het PHM stuk
   PHM hulp_ph(L,N);
   PHM G_b(L,N);

   //stop G(b) in G_b
   G_b.G(delta);

   //bereken G(rdm)^{-1}G(b)G(rdm)^{-1} en stop in hulp_ph
   hulp_ph.L_map(S.getG(),G_b);

   //tenslotte nog de antisymmetrische G hierop:
   TPM hulp(L,N);
   hulp.G(hulp_ph);

   //en optellen bij this
   *this += hulp;
#endif

   (*this) *= t;

   Proj_E(lineq);
}

double TPM::line_search(double t, SUP &S, const TPM &ham) const
{
   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   S.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta(L,N);

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp(L,N);

   hulp.L_map(S,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c = 0;

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance)
   {
      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;
   }

   return c;
}


/**
 * Reads a rdm from a HDF5 file stored as a full RDM:
 * no split up between a LxL block and a vector. This method
 * will read the full rdm and only keep the LxL block and the vector
 * @param filename the file to read
 */
void TPM::ReadFromFileFull(std::string filename)
{
   hid_t       file_id, dataset_id;
   herr_t      status;

   std::unique_ptr<Matrix> rdm_tmp(new Matrix (L*(2*L-1))); 

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   HDF5_STATUS_CHECK(file_id);

   dataset_id = H5Dopen(file_id, "/matrix", H5P_DEFAULT);
   HDF5_STATUS_CHECK(dataset_id);

   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdm_tmp->gMatrix());
   HDF5_STATUS_CHECK(status);

   status = H5Dclose(dataset_id);
   HDF5_STATUS_CHECK(status);

   status = H5Fclose(file_id);
   HDF5_STATUS_CHECK(status);

   for(int i=0;i<L;i++)
      for(int j=i;j<L;j++)
         getMatrix(0)(i,j) = getMatrix(0)(j,i) = (*rdm_tmp)(i,j);

   for(int i=0;i<gdimVector(0);i++)
      getVector(0)[i] = (*rdm_tmp)(L+i,L+i);
}

/*  
 * project the TPM object onto a space where (option == 0) Tr (Gamma' E) = 0, or (option == 1) Tr (Gamma' E) = e
 * @param lineq The object containing the linear constraints
 * @param option project onto (option = 0) 0 or (option = 1) e
 */
void TPM::Proj_E(const Lineq &lineq, int option)
{
   double brecht;

   for(int i=0;i<lineq.gnr();++i)
   {
      if(option == 1)
         //Tr(Gamma E~) - e~
         brecht = this->ddot(lineq.gE_ortho(i)) - lineq.ge_ortho(i);
      else
         brecht = this->ddot(lineq.gE_ortho(i));

      this->daxpy(-brecht,lineq.gE_ortho(i));
   }
}

/**
 * Created a vector with all necessary DOCI consistency constrains.
 * @return a vector of all constrains
 */
std::vector<TPM> TPM::DOCI_constrains() const
{
   std::vector<TPM> list;
   list.reserve(L);

   for(int a=0;a<L;a++)
   {
      TPM constr(L,N);
      constr = 0;

      constr(0,a,a) = -1;

      for(int b=0;b<L;b++)
         if(a!=b)
         {
            int i = (*s2t)(a,b);
            // divide by 4 to compensate for the degeneracy
            constr(0,i-L) = 1.0/(N/2.0-1)/4.0;
         }

      list.push_back(std::move(constr));
   }

   return list;
}

/**
 * Get an element from the diag vector part. Warning:
 * no check on a!=b. 
 * @param a first sp index
 * @param b second sp index
 * @return the element a b ; a b
 */
double TPM::getDiag(int a, int b) const
{
   assert(a<L && b<L);

   if(a==b)
      return 0;

   int i = (*s2t)(a,b);
   return (*this)(0,i-L);
}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param tpm_d the input TPM
 */
void TPM::S(const TPM &tpm_d, bool inverse)
{
   double a = 1.0;
   double b = 0.0;
   double c = 0.0;
   const int M = 2*L;

#ifdef __Q_CON

   a += 1.0;
   b += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   c += (2.0*N - M)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __G_CON

   a += 4.0;
   c += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   a += M - 4.0;
   b += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   c -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON
   
   a += 5.0*M - 8.0;
   b += 2.0/(N - 1.0);
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

   this->Q(a,b,c,tpm_d,inverse);
}

/**
 * Calculate the inverse Overlap matrix S^-1 a = b by
 * solving S b = a with conjugate gradient.
 * Store the result in *this.
 * @param b the input matrix
 * @param lineq the constrains to use
 * @returns the number of iterations
 */
int TPM::InverseS(TPM &b, const Lineq &lineq)
{
   *this = 0;

   //de r initialiseren op b
   TPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:
   TPM Hb(L,N);

   int cg_iter = 0;

   while(rr > 1.0e-10)
   {
      ++cg_iter;

      Hb.S(b);
      Hb.Proj_E(lineq);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;
   }

   return cg_iter;
}

/**
 * The down image of the DOCI-G image.
 * Fills the current object with the G down image of phm
 * @param phm the G matrix to use
 */
void TPM::G(const PHM &phm)
{
   std::vector<double> B11(L,0);
   std::vector<double> B22(L,0);

   for(int a=0;a<L;a++)
      for(int b=a+1;b<L;b++)
      {
         auto &block = phm.getBlock(a,b);

         B11[a] += block(0,0);
         B22[b] += block(1,1);
      }

   for(int i=0;i<gdimVector(0);++i)
   {
      const int a = (*t2s)(L+i,0);
      const int b = (*t2s)(L+i,1);

      auto &block = phm.getBlock(a,b);

      (*this)(0,i) = 0.25*(2.0/(N-1.0)*(phm[0](a,a)+phm[0](b,b)+B11[a]+B11[b]+B22[a]+B22[b]) - block(0,0) - block(1,1) + 2*phm[0](a,b));
   }

   for(int a=0;a<L;a++)
   {
      for(int b=a+1;b<L;b++)
      {
         auto &block = phm.getBlock(a,b);
         (*this)(0,a,b) = (*this)(0,b,a) = -block(0,1);
      }

      (*this)(0,a,a) = 1.0/(N-1.0)*(B11[a]+B22[a]+phm[0](a,a));
   }
}

/**
 * line search for interpolation
 * @param t barriere hight
 * @param rdm current rdm
 * @param ham the hamiltonian
 * @return the step size
 */
double TPM::line_search(double t, const TPM &rdm, const TPM &ham) const
{
   SUP P(L,N);

   P.fill(rdm);
   P.invert();

   return line_search(t,P,ham);
}

/**
 * Build pairing hamiltonian
 * @param g the pairing strength
 */
void TPM::pairing(double g)
{
   std::vector<double> Elevels(L);

   // picket fence sp levels
   for(int i=0;i<L;i++)
      Elevels[i] = -1.0/L*(L-i);

   std::vector<double> x(L);

   //pairing interaction term
   for(int a=0;a<L;a++)
      x[a] = 1.0;

   //normeren op 1/2
   double bright = 0.0;

   for(auto &elem: x)
      bright += elem*elem;

   bright *= 2;

   bright = 1.0/std::sqrt(bright);

   for(auto &elem: x)
      elem *= bright;

   // make our life easier
   auto calc_elem = [this,&g,&Elevels,&x] (int i, int j) {
      int a = (*t2s)(i,0);
      int b = (*t2s)(i,1);
      int c = (*t2s)(j,0);
      int d = (*t2s)(j,1);

      int a_ = a % L;
      int b_ = b % L;
      int c_ = c % L;
      int d_ = d % L;

      double result = 0;

      if(i==j)
      {
         result += (Elevels[a_] + Elevels[b_])*1.0/(N-1.0);

         if(a_ == b_)
            result -= 2 * g * x[a_] * x[a_];
      }
      else if(a_ == b_ && c_ == d_)
         result -= 2 * g * x[a_] * x[c_];

      return result;
   };

   auto& hamB = getMatrix(0);

   for(int i=0;i<L;++i)
      for(int j=i;j<L;++j)
         hamB(i,j) = hamB(j,i) = calc_elem(i,j);

   auto& hamV = getVector(0);

   for(int i=0;i<hamV.gn();i++)
      // keep in mind that the degen of the vector is 4. We need prefactor of 2, so
      // we end up with 0.5
      hamV[i] = 0.5*calc_elem(L+i,L+i) + 0.5*calc_elem(L*L+i,L*L+i);
}

/**
 * Rotate this TPM object with a jacobi rotation between orbitals
 * k and l over an angle of angle in the DOCI space
 * @param k the first orbital
 * @param l the second orbital
 * @param angle the angle to rotate over
 */
void TPM::rotate_doci(int k, int l, double angle)
{
   TPM new_rdm(*this);

   auto& rdmB = new_rdm.getMatrix(0);
   auto& rdmV = new_rdm.getVector(0);

   const double cos = std::cos(angle);
   const double sin = std::sin(angle);
   const double cos2 = cos*cos;
   const double sin2 = sin*sin;
   const double cos4 = cos2*cos2;
   const double sin4 = sin2*sin2;
   const double cos2sin2 = cos2*sin2;

   for(int p=0;p<L;p++)
      if(p == k || p ==l)
         continue;
      else
      {
         rdmB(k,p) = rdmB(p,k) = cos2*(*this)(0,k,p)+sin2*(*this)(0,l,p);

         rdmB(l,p) = rdmB(p,l) = cos2*(*this)(0,l,p)+sin2*(*this)(0,k,p);

         int idx1 = (*s2t)(k,p) - L;
         int idx2 = (*s2t)(l,p) - L;

         rdmV[idx1] = cos2 * (*this)(k,p,k,p) + sin2 * (*this)(l,p,l,p);
         rdmV[idx2] = cos2 * (*this)(l,p,l,p) + sin2 * (*this)(k,p,k,p);
      }

   // k \bar k; k \bar k
   rdmB(k,k) = cos4 * (*this)(0,k,k) + 2 * cos2sin2 * (*this)(0,k,l) + sin4 * (*this)(0,l,l) + \
               2 * cos2sin2 * (*this)(k,l,k,l);

   // l \bar l; l \bar l
   rdmB(l,l) = cos4 * (*this)(0,l,l) + 2 * cos2sin2 * (*this)(0,k,l) + sin4 * (*this)(0,k,k) + \
               2 * cos2sin2 * (*this)(k,l,k,l);

   // k \bar k; l \bar l
   rdmB(k,l) = cos2sin2 * (*this)(0,k,k) + (cos4 + sin4) * (*this)(0,k,l) + cos2sin2 * (*this)(0,l,l) - \
               2 * cos2sin2 * (*this)(k,l,k,l);
   rdmB(l,k) = rdmB(k,l);

   int idx = (*s2t)(k,l) - L;
   rdmV[idx] += cos2sin2*((*this)(0,k,k)+(*this)(0,l,l)-2*(*this)(0,k,l))+(cos4+sin4)*(*this)(k,l,k,l);
   rdmV[idx] *= 0.5;

   (*this) = std::move(new_rdm);
}

/**
 * Calculate the energy change when you rotate orbital k and l over an angle of theta
 * in the DOCI space.
 * @param ham the hamiltonian to use
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotate over
 * @return the new energy
 */
double TPM::calc_rotate_doci(const TPM &ham, int k, int l, double theta) const
{
   double energy = 0;

   auto A = [&] (int a, int b, int a2, int b2) -> double { return ham(0,a,b)*(*this)(0,a2,b2) + 2*ham(a,b,a,b)*(*this)(a2,b2,a2,b2); };
   auto B = [&] (int a, int b, int a2, int b2) -> double { return ham(0,a,b)*(*this)(a2,b2,a2,b2) + ham(a,b,a,b)*((*this)(0,a2,b2) - (*this)(a2,b2,a2,b2)); };

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      for(int b=0;b<L;b++)
         if(b==k || b==l)
            continue;
         else
            energy += A(a,b,a,b);

      energy += 2 * std::cos(theta)*std::cos(theta) * ( A(k,a,k,a) + A(l,a,l,a) ) + 2 * std::sin(theta)*std::sin(theta) * ( A(k,a,l,a) + A(l,a,k,a) );
   }

   // cos^4
   energy += std::cos(theta)*std::cos(theta)*std::cos(theta)*std::cos(theta) * (A(l,l,l,l)+A(k,k,k,k)+2*A(k,l,k,l));

   // sin^4
   energy += std::sin(theta)*std::sin(theta)*std::sin(theta)*std::sin(theta) * (A(k,k,l,l)+A(l,l,k,k)+2*A(k,l,k,l));

   // sin^2 cos^2
   energy += std::cos(theta)*std::cos(theta)*std::sin(theta)*std::sin(theta) * 2.0 * (A(l,l,k,l)+A(k,l,l,l)+A(k,l,k,k)+A(k,k,k,l) + \
         B(l,l,k,l)+B(k,l,l,l)+B(k,l,k,k)+B(k,k,k,l) -2*B(k,l,k,l));

   return energy;
}

/**
 * The slow version of TPM::calc_rotate_doci()
 * @param ham the hamiltonian to use
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotate over
 * @return the new energy
 */
double TPM::calc_rotate_slow_doci(const TPM &ham, int k, int l, double theta) const
{
   double energy = 0;

   Matrix rot(2*L);
   rot = 0;
   for(int i=0;i<2*L;i++)
      rot(i,i) = 1.0;
   rot(k,k) = std::cos(theta);
   rot(l,l) = std::cos(theta);
   rot(k,l) = -1.0*std::sin(theta);
   rot(l,k) = std::sin(theta);

   rot(k+L,k+L) = std::cos(theta);
   rot(l+L,l+L) = std::cos(theta);
   rot(k+L,l+L) = -1.0*std::sin(theta);
   rot(l+L,k+L) = std::sin(theta);

   for(int a=0;a<2*L;a++)
      for(int b=0;b<2*L;b++)
         for(int c=0;c<2*L;c++)
            for(int d=0;d<2*L;d++)
            {
               if(fabs((*this)(a,b,c,d)) < 1e-14)
                  continue;

               for(int a2=0;a2<2*L;a2++)
                  for(int b2=0;b2<2*L;b2++)
                     for(int c2=0;c2<2*L;c2++)
                        for(int d2=0;d2<2*L;d2++)
                        {
                           if(fabs(ham(a2,b2,c2,d2)) < 1e-14)
                              continue;

                           energy += 0.25 * ham(a,b,c,d) * rot(a,a2) * rot(b,b2) * \
                                     rot(c,c2)* rot(d,d2)* (*this)(a2,b2,c2,d2);
                        }
            }

   return energy;
}


/**
 * Find the minimum with Newton-Raphson for the angle of a jacobi rotation between
 * orbitals k and l in the DOCI space.
 * @param ham the hamiltonian to use
 * @param k the first orbital
 * @param l the second orbital
 * @param start_angle the starting point for the Newton-Raphson (defaults to zero)
 * @return pair of the angle with the lowest energy and boolean, true => minimum, false => maximum
 */
std::pair<double,bool> TPM::find_min_angle_doci(const TPM &ham, int k, int l, double start_angle) const
{
   double theta = start_angle;

   auto A = [&] (int a, int b, int a2, int b2) -> double { return ham(a,a+L,b,b+L)*(*this)(a2,a2+L,b2,b2+L) + 2*ham(a,b,a,b)*(*this)(a2,b2,a2,b2); };
   auto B = [&] (int a, int b, int a2, int b2) -> double { return ham(a,a+L,b,b+L)*(*this)(a2,b2,a2,b2) + ham(a,b,a,b)*((*this)(a2,a2+L,b2,b2+L) - (*this)(a2,b2,a2,b2)); };


   double sum1 = 0;
   double sum2 = 0;

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      sum1 += A(k,a,k,a) + A(l,a,l,a);
      sum2 += A(k,a,l,a) + A(l,a,k,a);
   }

   const double cos4 = A(l,l,l,l) + A(k,k,k,k) + 2*A(k,l,k,l);
   const double sin4 = A(k,k,l,l) + A(l,l,k,k) + 2*A(k,l,k,l);
   const double cos2sin2 = A(l,l,k,l) + A(k,l,l,l) + A(k,l,k,k) + A(k,k,k,l) + \
                           B(l,l,k,l) + B(k,l,l,l) + B(k,l,k,k) + B(k,k,k,l) - 2*B(l,k,l,k);

   auto gradient = [&cos4,&sin4,&cos2sin2,&sum1,&sum2] (double t) -> double {
      return -4*std::sin(t)*std::cos(t) * (std::cos(t)*std::cos(t)*(cos4 + sin4 - 2 * cos2sin2) + cos2sin2 - sin4 + sum1 - sum2);
   };

   auto hessian = [&cos4,&sin4,&cos2sin2,&sum1,&sum2] (double t) -> double {
      return ((-cos4-sin4+2*cos2sin2) * 16 * std::cos(t)*std::cos(t) + (12*cos4+20*sin4-8*sum1+8*sum2-32*cos2sin2))*std::cos(t)*std::cos(t)-4*sin4+4*sum1-4*sum2+4*cos2sin2;
   };

   int max_iters = 100;
   double new_theta;

   for(int iter=0;iter<max_iters;iter++)
   {
      new_theta = theta - gradient(theta)/hessian(theta);

      theta = new_theta;

      if(fabs(gradient(theta)) < 1e-12)
         break;
   }

//   int Na = 5000;
//   for(int i=0;i<Na;i++)
//   {
//      double t = 2.0 * M_PI / (1.0*Na) * i;
//      std::cout << t << "\t" << calc_rotate(ham,k,l,t)  << "\t" << gradient(t) << "\t" << hessian(t) << std::endl;
//// calc_rotate_slow(ham,k,l,t) << "\t" << 
//   }

   return std::make_pair(theta, hessian(theta)>0);
}

/**
 * The slow version of TPM::calc_rotate()
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotate over
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 * @return the new energy
 */
double TPM::calc_rotate_slow(int k, int l, double theta, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const
{
   assert(k!=l);

   double energy = 0;

   Matrix rot(L);
   rot = 0;
   rot.unit();
   rot(k,k) = std::cos(theta);
   rot(l,l) = std::cos(theta);
   rot(k,l) = -1*std::sin(theta);
   rot(l,k) = std::sin(theta);

   for(int a=0;a<L;a++)
      for(int b=0;b<L;b++)
         for(int a2=0;a2<L;a2++)
            for(int b2=0;b2<L;b2++)
            {
               energy += 2.0/(N-1.0) * (rot(a,a2)*rot(a,b2)+rot(b,a2)*rot(b,b2))*(*this)(a,b,a,b)*T(a2,b2);

               if(a==b)
                  energy += 2.0/(N-1.0)*rot(a,a2)*rot(a,b2)*(*this)(a,a+L,a,a+L)*T(a2,b2);

               for(int c2=0;c2<L;c2++)
                  for(int d2=0;d2<L;d2++)
                     energy += V(a2,b2,c2,d2) * ( \
                           rot(a,a2)*rot(a,b2)*rot(b,c2)*rot(b,d2)* (*this)(a,a+L,b,b+L) + \
                           rot(a,a2)*rot(b,b2)*( \
                              2*rot(a,c2)*rot(b,d2) -
                              rot(b,c2)*rot(a,d2) ) * \
                           (*this)(a,b,a,b) \
                           );
            }

   return energy;
}

/**
 * Calculate the energy change when you rotate orbital k and l over an angle of theta with
 * the rotation in the full space of the orbitals.
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotate over
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 * @return the new energy
 */
double TPM::calc_rotate(int k, int l, double theta, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const
{
   assert(k!=l);

   const TPM &rdm = *this;

   double energy = 4/(N-1.0)*(T(k,k)+T(l,l)) * rdm(k,l,k,l);

   double cos2 = 2.0/(N-1.0)*(T(k,k)*rdm(k,k+L,k,k+L)+T(l,l)*rdm(l,l+L,l,l+L));

   double sin2 = 2.0/(N-1.0)*(T(l,l)*rdm(k,k+L,k,k+L)+T(k,k)*rdm(l,l+L,l,l+L));

   // 2sincos actually
   double sincos = 2.0/(N-1.0)*T(k,l)*(rdm(l,l+L,l,l+L)-rdm(k,k+L,k,k+L));

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      energy += 2.0/(N-1.0) * T(a,a) * (rdm(a,a+L,a,a+L)+2*rdm(a,k,a,k)+2*rdm(a,l,a,l));

      for(int b=0;b<L;b++)
      {
         if(b==k || b==l)
            continue;

         energy += 2.0/(N-1.0) * (T(a,a)+T(b,b)) * rdm(a,b,a,b);

         energy += V(a,a,b,b) * rdm(a,a+L,b,b+L);

         energy += (2*V(a,b,a,b)-V(a,b,b,a)) * rdm(a,b,a,b);
      }

      cos2 += 2*V(k,k,a,a)*rdm(k,k+L,a,a+L)+2*V(l,l,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(k,a,k,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(l,a,l,a);

      sin2 += 2*V(l,l,a,a)*rdm(k,k+L,a,a+L)+2*V(k,k,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(l,a,l,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(k,a,k,a);

      sincos += 2*V(k,l,a,a)*(rdm(l,l+L,a,a+L)-rdm(k,k+L,a,a+L))+2*(2*V(k,a,l,a)-V(k,a,a,l)+2.0/(N-1.0)*T(k,l))*(rdm(l,a,l,a)-rdm(k,a,k,a));
   }

   const double cos4 = V(k,k,k,k)*rdm(k,k+L,k,k+L)+V(l,l,l,l)*rdm(l,l+L,l,l+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   const double sin4 = V(k,k,k,k)*rdm(l,l+L,l,l+L)+V(l,l,l,l)*rdm(k,k+L,k,k+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   // 2 x
   const double cos2sin2 = (2*V(k,k,l,l)+V(k,l,k,l))*(rdm(k,k+L,k,k+L)+rdm(l,l+L,l,l+L))+((V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l))))*rdm(k,k+L,l,l+L)+(V(k,k,k,k)+V(l,l,l,l)-6*V(k,k,l,l)+2*V(k,l,k,l))*rdm(k,l,k,l);

   // 4 x
   const double sin3cos = V(k,l,k,k)*rdm(l,l+L,l,l+L)-V(k,l,l,l)*rdm(k,k+L,k,k+L)-(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // 4 x
   const double cos3sin = V(k,l,l,l)*rdm(l,l+L,l,l+L)-V(k,l,k,k)*rdm(k,k+L,k,k+L)+(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   const double cos = std::cos(theta);
   const double sin = std::sin(theta);

   energy += cos*cos*cos*cos*cos4;

   energy += sin*sin*sin*sin*sin4;

   energy += cos*cos*cos2;

   energy += sin*sin*sin2;

   energy += 2*sin*cos*sincos;

   energy += 2*cos*cos*sin*sin*cos2sin2;

   energy += 4*cos*sin*sin*sin*sin3cos;

   energy += 4*cos*cos*cos*sin*cos3sin;




//   const double cos = std::cos(theta);
//   const double sin = std::sin(theta);
//   const double cos2 = cos*cos;
//   const double sin2 = sin*sin;
//   const double cos4 = cos2*cos2;
//   const double sin4 = sin2*sin2;
//   const double cossin = cos*sin;
//   const double cos2sin2 = cos2*sin2;
//   const double cos3sin = cos2*cossin;
//   const double cossin3 = cossin*sin2;
//
//   const double fac1 = cos2*T(k,k)+sin2*T(l,l)-2*cossin*T(k,l);
//   const double fac2 = cos2*T(l,l)+sin2*T(k,k)+2*cossin*T(k,l);
//
//   for(int a=0;a<L;a++)
//   {
//      if(a==k || a==l)
//         continue;
//
//      energy += 2.0/(N-1.0) * T(a,a) * (*this)(a,a+L,a,a+L);
//
//      for(int b=0;b<L;b++)
//      {
//         if(b==k || b==l)
//            continue;
//
//         energy += 2.0/(N-1.0) * (T(a,a)+T(b,b)) * (*this)(a,b,a,b);
//
//         energy += V(a,a,b,b) * (*this)(a,a+L,b,b+L);
//
//         energy += (2*V(a,b,a,b)-V(a,b,b,a)) * (*this)(a,b,a,b);
//      }
//
//      energy += 4.0/(N-1.0) * (T(a,a)+ fac1) * (*this)(a,k,a,k);
//      energy += 4.0/(N-1.0) * (T(a,a)+ fac2) * (*this)(a,l,a,l);
//
//      energy += 2*(cos2*V(k,k,a,a)-2*cossin*V(k,l,a,a)+sin2*V(l,l,a,a))*(*this)(k,k+L,a,a+L);
//
//      energy += 2*(cos2*V(l,l,a,a)+2*cossin*V(k,l,a,a)+sin2*V(k,k,a,a))*(*this)(l,l+L,a,a+L);
//
//      energy += 2*(cos2*(2*V(k,a,k,a)-V(k,a,a,k))-2*cossin*(2*V(k,a,l,a)-V(k,a,a,l))+sin2*(2*V(l,a,l,a)-V(l,a,a,l)))*(*this)(k,a,k,a);
//
//      energy += 2*(cos2*(2*V(l,a,l,a)-V(l,a,a,l))+2*cossin*(2*V(l,a,k,a)-V(k,a,a,l))+sin2*(2*V(k,a,k,a)-V(k,a,a,k)))*(*this)(l,a,l,a);
//   }
//
//   energy += 2.0/(N-1.0)*fac1*(*this)(k,k+L,k,k+L);
//   energy += 2.0/(N-1.0)*fac2*(*this)(l,l+L,l,l+L);
//   energy += 4.0/(N-1.0)*(T(k,k)+T(l,l))*(*this)(k,l,k,l);
//
//   energy += 2*(cos2sin2*(V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l)))+(cos4+sin4)*V(k,k,l,l)+2*(cos3sin-cossin3)*(V(k,l,k,k)-V(k,l,l,l)))*(*this)(k,k+L,l,l+L);
//
//   energy += (cos4*V(k,k,k,k)+sin4*V(l,l,l,l)+cos2sin2*(4*V(k,k,l,l)+2*V(k,l,k,l))-4*cossin3*V(k,l,l,l)-4*cos3sin*V(k,l,k,k))*(*this)(k,k+L,k,k+L);
//
//   energy += (sin4*V(k,k,k,k)+cos4*V(l,l,l,l)+cos2sin2*(4*V(k,k,l,l)+2*V(k,l,k,l))+4*cossin3*V(k,l,k,k)+4*cos3sin*V(k,l,l,l))*(*this)(l,l+L,l,l+L);
//
//   energy += 2*(cos2sin2*(V(k,k,k,k)+V(l,l,l,l)-6*V(k,k,l,l)+2*V(k,l,k,l))+(cos4+sin4)*(2*V(k,l,k,l)-V(k,k,l,l))+2*(cos3sin-cossin3)*(V(k,l,k,k)-V(k,l,l,l)))*(*this)(k,l,k,l);

   return energy;
}

/**
 * Find the minimum with Newton-Raphson for the angle of a jacobi rotation between
 * orbitals k and l in the DOCI space.
 * @param k the first orbital
 * @param l the second orbital
 * @param start_angle the starting point for the Newton-Raphson (defaults to zero)
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 * @return pair of the angle with the lowest energy and boolean, true => minimum, false => maximum
 */
std::pair<double,bool> TPM::find_min_angle(int k, int l, double start_angle, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const
{
   assert(k!=l);

   double theta = start_angle;

   const TPM &rdm = *this;

   double cos2 = 2.0/(N-1.0)*(T(k,k)*rdm(k,k+L,k,k+L)+T(l,l)*rdm(l,l+L,l,l+L));

   double sin2 = 2.0/(N-1.0)*(T(l,l)*rdm(k,k+L,k,k+L)+T(k,k)*rdm(l,l+L,l,l+L));

   // 2sincos actually
   double sincos = 2.0/(N-1.0)*T(k,l)*(rdm(l,l+L,l,l+L)-rdm(k,k+L,k,k+L));

   for(int a=0;a<L;a++)
   {
      if(a==k || a==l)
         continue;

      cos2 += 2*V(k,k,a,a)*rdm(k,k+L,a,a+L)+2*V(l,l,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(k,a,k,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(l,a,l,a);

      sin2 += 2*V(l,l,a,a)*rdm(k,k+L,a,a+L)+2*V(k,k,a,a)*rdm(l,l+L,a,a+L)+2*(2*V(k,a,k,a)-V(k,a,a,k)+2.0/(N-1.0)*T(k,k))*rdm(l,a,l,a)+2*(2*V(l,a,l,a)-V(l,a,a,l)+2.0/(N-1.0)*T(l,l))*rdm(k,a,k,a);

      sincos += 2*V(k,l,a,a)*(rdm(l,l+L,a,a+L)-rdm(k,k+L,a,a+L))+2*(2*V(k,a,l,a)-V(k,a,a,l)+2.0/(N-1.0)*T(k,l))*(rdm(l,a,l,a)-rdm(k,a,k,a));
   }

   const double cos4 = V(k,k,k,k)*rdm(k,k+L,k,k+L)+V(l,l,l,l)*rdm(l,l+L,l,l+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   const double sin4 = V(k,k,k,k)*rdm(l,l+L,l,l+L)+V(l,l,l,l)*rdm(k,k+L,k,k+L)+2*V(k,k,l,l)*rdm(k,k+L,l,l+L)+2*(2*V(k,l,k,l)-V(k,k,l,l))*rdm(k,l,k,l);

   // 2 x
   const double cos2sin2 = (2*V(k,k,l,l)+V(k,l,k,l))*(rdm(k,k+L,k,k+L)+rdm(l,l+L,l,l+L))+((V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l))))*rdm(k,k+L,l,l+L)+(V(k,k,k,k)+V(l,l,l,l)-6*V(k,k,l,l)+2*V(k,l,k,l))*rdm(k,l,k,l);

   // 4 x
   const double sin3cos = V(k,l,k,k)*rdm(l,l+L,l,l+L)-V(k,l,l,l)*rdm(k,k+L,k,k+L)-(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // 4 x
   const double cos3sin = V(k,l,l,l)*rdm(l,l+L,l,l+L)-V(k,l,k,k)*rdm(k,k+L,k,k+L)+(V(k,l,k,k)-V(k,l,l,l))*(rdm(k,k+L,l,l+L)+rdm(k,l,k,l));

   // A*cos(t)^4+B*sin(t)^4+C*cos(t)^2+D*sin(t)^2+2*E*cos(t)*sin(t)+2*F*cos(t)^2*sin(t)^2+4*G*sin(t)*cos(t)^3+4*H*sin(t)^3*cos(t)

   // (16*G-16*H)*cos(t)^4+(-4*A-4*B+8*F)*sin(t)*cos(t)^3+(4*E-12*G+20*H)*cos(t)^2+(4*B-2*C+2*D-4*F)*sin(t)*cos(t)-2*E-4*H
   auto gradient = [&] (double theta) -> double { 
      double cos = std::cos(theta);
      double sin = std::sin(theta);

      double res = 16*(cos3sin-sin3cos)*cos*cos*cos*cos - 4*(cos4+sin4-2*cos2sin2)*sin*cos*cos*cos + 4*(sincos-3*cos3sin+5*sin3cos)*cos*cos + 2*(2*sin4-cos2+sin2-2*cos2sin2)*sin*cos-2*sincos-4*sin3cos;

      return res;
   };

   // (-16*A-16*B+32*F)*cos(t)^4+(-64*G+64*H)*sin(t)*cos(t)^3+(12*A+20*B-4*C+4*D-32*F)*cos(t)^2+(-8*E+24*G-40*H)*sin(t)*cos(t)-4*B+2*C-2*D+4*F
   auto hessian = [&] (double theta) -> double { 
      double cos = std::cos(theta);
      double sin = std::sin(theta);

      double res = -16*(cos4+sin4-2*cos2sin2)*cos*cos*cos*cos+64*(sin3cos-cos3sin)*sin*cos*cos*cos+4*(3*cos4+5*sin4-cos2+sin2-8*cos2sin2)*cos*cos-8*(sincos-3*cos3sin+5*sin3cos)*sin*cos+2*(cos2-2*sin4-sin2+2*cos2sin2);
      
      return res;
   };

//   int Na = 5000;
//   for(int i=0;i<Na;i++)
//   {
//      double t = 2.0 * M_PI / (1.0*Na) * i;
//      std::cout << t << "\t" << calc_rotate(k,l,t,T,V)  << "\t" << gradient(t) << "\t" << hessian(t) << std::endl;
//   }

   const int max_iters = 20;
   const double convergence = 1e-12;

   double change = gradient(theta)*theta+hessian(theta)*theta*theta/2.0;

   // if it goes uphill, try flipping sign and try again
   if(change>0)
      theta *= -1;

   for(int iter=0;iter<max_iters;iter++)
   {
      double dx = gradient(theta)/hessian(theta);

      theta -= dx;

      if(fabs(dx) < convergence)
         break;
   }

   if(hessian(theta)<0)
      std::cout << "Found max!" << std::endl;

   return std::make_pair(theta, hessian(theta)>0);
}

/**
 * Rotate this TPM object with a jacobi rotation between orbitals
 * k and l over angle in the full space.
 * @param k the first orbital
 * @param l the second orbital
 * @param angle the angle to rotate over
 * @param T function that returns the one-particle matrix elements
 * @param V function that returns the two-particle matrix elements
 */
void TPM::rotate(int k, int l, double angle, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V)
{
   assert(k!=l);

   auto& rdmB = getMatrix(0);
   auto& rdmV = getVector(0);

   const double cos = std::cos(angle);
   const double sin = std::sin(angle);
   const double cos2 = cos*cos;
   const double sin2 = sin*sin;
   const double cos4 = cos2*cos2;
   const double sin4 = sin2*sin2;
   const double cossin = cos*sin;
   const double cos2sin2 = cos2*sin2;
   const double cos3sin = cos2*cossin;
   const double cossin3 = cossin*sin2;

   for(int p=0;p<L;p++)
   {
      if(p == k || p ==l)
         continue;

      // k \bar k ; p \bar p
      rdmB(k,p) = rdmB(p,k) = cos2*V(k,k,p,p)-2*cossin*V(k,l,p,p)+sin2*V(l,l,p,p);
      // l \bar l ; p \bar p
      rdmB(l,p) = rdmB(p,l) = cos2*V(l,l,p,p)+2*cossin*V(k,l,p,p)+sin2*V(k,k,p,p);

      int idx = (*s2t)(k,p) - L;

      // k p ; k p
      rdmV[idx] = 1.0/(N-1.0) * (T(p,p) + cos2*T(k,k)-2*cossin*T(k,l)+sin2*T(l,l)); 
      rdmV[idx] += cos2*(V(k,p,k,p)-0.5*V(k,p,p,k))-2*cossin*(V(k,p,l,p)-0.5*V(k,p,p,l))+sin2*(V(l,p,l,p)-0.5*V(l,p,p,l));

      idx = (*s2t)(l,p) - L;

      // l p ; l p
      rdmV[idx] = 1.0/(N-1.0) * (T(p,p) + cos2*T(l,l)+2*cossin*T(k,l)+sin2*T(k,k)); 
      rdmV[idx] += cos2*(V(l,p,l,p)-0.5*V(l,p,p,l))+2*cossin*(V(l,p,k,p)-0.5*V(k,p,p,l))+sin2*(V(k,p,k,p)-0.5*V(k,p,p,k));
   }

   // k \bar k ; k \bar k
   rdmB(k,k) = 2.0/(N-1.0) * (cos2*T(k,k)-2*cossin*T(k,l)+sin2*T(l,l));
   rdmB(k,k) += cos4*V(k,k,k,k)+sin4*V(l,l,l,l)+cos2sin2*(4*V(k,k,l,l)+2*V(k,l,k,l))-4*cossin3*V(k,l,l,l)-4*cos3sin*V(k,l,k,k);

   // l \bar l ; l \bar l
   rdmB(l,l) = 2.0/(N-1.0) * (cos2*T(l,l)+2*cossin*T(k,l)+sin2*T(k,k));
   rdmB(l,l) += sin4*V(k,k,k,k)+cos4*V(l,l,l,l)+cos2sin2*(4*V(k,k,l,l)+2*V(k,l,k,l))+4*cossin3*V(k,l,k,k)+4*cos3sin*V(k,l,l,l);

   // k \bar k ; l \bar l
   rdmB(k,l) = rdmB(l,k) = cos2sin2*(V(k,k,k,k)+V(l,l,l,l)-2*(V(k,l,k,l)+V(k,k,l,l)))+(cos4+sin4)*V(k,k,l,l)+2*(cos3sin-cossin3)*(V(k,l,k,k)-V(k,l,l,l));

   // k l ; k l
   int idx = (*s2t)(k,l) - L;
   rdmV[idx] = 1.0/(N-1.0)*(T(k,k)+T(l,l)) + cos2sin2*(0.5*(V(k,k,k,k)+V(l,l,l,l))-3*V(k,k,l,l)+V(k,l,k,l))+(cos4+sin4)*(V(k,l,k,l)-0.5*V(k,k,l,l))+(cos3sin-cossin3)*(V(k,l,k,k)-V(k,l,l,l));
}

/*  vim: set ts=3 sw=3 expandtab :*/
