#include <assert.h>

#include "EIG.h"


EIG::EIG(BlockMatrix &mat): BlockVector(mat.gnr()) 
{
   for(int i=0;i<mat.gnr();i++)
   {
      setDim(i,mat.gdim(i),mat.gdeg(i));

      (*this)[i].diagonalize(mat[i]);
   }
}

EIG::EIG(Container &cont): BlockVector(cont.gnr())
{
   for(int i=0;i<cont.gnMatrix();i++)
   {
      setDim(i,cont.gdimMatrix(i),cont.gdegMatrix(i));

      (*this)[i].diagonalize(cont.getMatrix(i));
   }

   for(int i=0;i<cont.gnVector();i++)
   {
      setDim(i+cont.gnMatrix(), cont.gdimVector(i), cont.gdegVector(i));
      (*this)[i+cont.gnMatrix()] = cont.getVector(i);
   }
}

EIG::EIG(SUP &sup): BlockVector(sup.gnr())
{
   int tel = 0;

   for(int i=0;i<sup.getI().gnMatrix();i++)
   {
      setDim(tel,sup.getI().gdimMatrix(i),sup.getI().gdegMatrix(i));
      (*this)[tel++].diagonalize(sup.getI().getMatrix(i));
   }

   for(int i=0;i<sup.getI().gnVector();i++)
   {
      setDim(tel, sup.getI().gdimVector(i), sup.getI().gdegVector(i));
      (*this)[tel++] = sup.getI().getVector(i);
   }

#ifdef __Q_CON
   for(int i=0;i<sup.getQ().gnMatrix();i++)
   {
      setDim(tel,sup.getQ().gdimMatrix(i),sup.getQ().gdegMatrix(i));
      (*this)[tel++].diagonalize(sup.getQ().getMatrix(i));
   }

   for(int i=0;i<sup.getQ().gnVector();i++)
   {
      setDim(tel, sup.getQ().gdimVector(i), sup.getQ().gdegVector(i));
      (*this)[tel++] = sup.getQ().getVector(i);
   }
#endif

#ifdef __G_CON
   for(int i=0;i<sup.getG().gnr();i++)
   {
      setDim(tel, sup.getG().gdim(i), sup.getG().gdeg(i));
      (*this)[tel++].diagonalize(sup.getG()[i]);
   }
#endif
}

double EIG::min() const
{
   double min = (*this)[0].min();

   for(int i=1;i<gnr();i++)
      if(min > (*this)[i].min())
         min = (*this)[i].min();

   return min;
}

double EIG::max() const
{
   double max = (*this)[0].max();

   for(int i=1;i<gnr();i++)
      if(max < (*this)[i].max())
         max = (*this)[i].max();

   return max;
}

double EIG::lsfunc(double a) const
{
   double res = 0;

   for(int i=0;i<gnr();i++)
      for(int j=0;j<gdim(i);j++)
        res += gdeg(i) * (*this)(i,j)/(1+a*(*this)(i,j));

   return res;
}

/*  vim: set ts=3 sw=3 expandtab :*/
