#include <sstream>
#include <assert.h>
#include <hdf5.h>

#include "include.h"

// default empty
std::unique_ptr<helpers::matrix> TPM::s2t = nullptr;
std::unique_ptr<helpers::matrix> TPM::t2s = nullptr;

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

   s2t.reset(new helpers::matrix(M,M));
   (*s2t) = -1; // if you use something you shouldn't, this will case havoc

   t2s.reset(new helpers::matrix(n_tp,2));
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

std::ostream &operator<<(std::ostream &output,TPM &tpm)
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


   // make our life easier
   auto calc_elem = [this,&OEI,&TEI] (int i, int j) {
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
         result += (OEI(a_,a_) + OEI(b_,b_))/(N - 1.0);

      // tp terms:

      // a \bar a ; b \bar b
      if(b==(a+L) && d==(c+L))
         result += TEI(a_*L + b_, c_*L + d_);

      // a b ; a b
      if(i==j && a<L && b<L && a!=b)
         result += TEI(a_*L + b_, c_*L + d_) - TEI(a_*L + b_, d_*L + c_);

      // \bar a \bar b ; \bar a \bar b
      if(i==j && a>=L && b>=L && a!=b)
         result += TEI(a_*L + b_, c_*L + d_) - TEI(a_*L + b_, d_*L + c_);

      // a \bar b ; a \bar b
      if(i==j && a<L && b>=L && a%L!=b%L)
         result += TEI(a_*L + b_, c_*L + d_);

      // \bar a b ; \bar a b
      if(i==j && a>=L && b<L && a%L!=b%L)
         result += TEI(a_*L + b_, c_*L + d_);

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
   double spin = 0.0;

   for(int i = 0;i<n;++i)
   {
      int a = (*t2s)(i,0);
      int b = (*t2s)(i,1);

      int s_a = ( 1.0 - 2 * (a / L) )/2;
      int s_b = ( 1.0 - 2 * (b / L) )/2;

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

int TPM::solve(double t, const SUP &S, TPM &grad, const Lineq &lineq)
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

void TPM::G(const PHM &phm)
{
   SPM spm(L,N);
   spm.bar(1.0/(N-1.0), phm);

   for(int i=0;i<L;i++)
   {
      int a = (*t2s)(i,0);

      for(int j=i;j<L;j++)
      {
         int c = (*t2s)(j,0);

         (*this)(0,i,j) = 2*(phm(a+L,c+L,c,a) - phm(a,c+L,c,a+L));

         (*this)(0,j,i) = (*this)(0,i,j);
      }

      (*this)(0,i,i) += 2 * spm(0,a);
   }

   for(int i=0;i<gdimVector(0);i++)
   {
      int a = (*t2s)(L+i,0);
      int b = (*t2s)(L+i,1);

      (*this)(0,i) = spm(0,a) + spm(0,b);
      (*this)(0,i) += 2*phm(a,a,b,b) - phm(a,b,a,b) - phm(b,a,b,a);
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

/*  vim: set ts=3 sw=3 expandtab :*/
