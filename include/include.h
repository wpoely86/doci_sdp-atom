#ifndef INCLUDE_H
#define INCLUDE_H

#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif


#include "Matrix.h"
#include "Vector.h"
#include "BlockStructure.h"
#include "Container.h"
#include "TPM.h"

#endif

/*  vim: set ts=3 sw=3 expandtab :*/ 
