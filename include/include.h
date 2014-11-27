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
#include "Tools.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"

#include "SUP.h"
#include "EIG.h"

#include "Lineq.h"

#include "Method.h"

#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;


#endif

/*  vim: set ts=3 sw=3 expandtab :*/ 
