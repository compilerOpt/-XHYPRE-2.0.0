/*BHEADER**********************************************************************
 * 
 * XHYPRE 2020.10.29
 *
 * autor: Chuanying Li
 *
 ***********************************************************************EHEADER*/
 
#ifndef hypre_HIG_MV_HEADER
#define hypre_HIG_MV_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HYPRE_config.h>

#include "seq_mv.h"
#include "_hypre_struct_mv.h"
//#include "protos.h"

#include "_hypre_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif


/* hig_csr_matvec.c */
// y[offset:end] = alpha*A[offset:end,:]*x + beta*b[offset:end]
HYPRE_Int hypre_HIGCSRMatrixMatvecOutOfPlace ( HYPRE_Complex alpha , hypre_CSRMatrix *A , hypre_Vector *x , HYPRE_Complex beta , hypre_Vector *b, hypre_Vector *y, HYPRE_Int offset );
// y = alpha*A + beta*y
HYPRE_Int hypre_HIGCSRMatrixMatvec ( HYPRE_Complex alpha , hypre_CSRMatrix *A , hypre_Vector *x , HYPRE_Complex beta , hypre_Vector *y );
HYPRE_Int hypre_HIGCSRMatrixMatvecT ( HYPRE_Complex alpha , hypre_CSRMatrix *A , hypre_Vector *x , HYPRE_Complex beta , hypre_Vector *y );
HYPRE_Int hypre_HIGCSRMatrixMatvec_FF ( HYPRE_Complex alpha , hypre_CSRMatrix *A , hypre_Vector *x , HYPRE_Complex beta , hypre_Vector *y , HYPRE_Int *CF_marker_x , HYPRE_Int *CF_marker_y , HYPRE_Int fpt );

/* hig_vector.c */
HYPRE_Real hypre_HIGSeqVectorInnerProd ( hypre_Vector *x , hypre_Vector *y );

/* struct_hig_matvec.c */
HYPRE_Int hypre_StructHIGMatvecCompute ( void *matvec_vdata , HYPRE_Complex alpha , hypre_StructMatrix *A , hypre_StructVector *x , HYPRE_Complex beta , hypre_StructVector *y );
HYPRE_Int hypre_StructHIGMatvec ( HYPRE_Complex alpha , hypre_StructMatrix *A , hypre_StructVector *x , HYPRE_Complex beta , hypre_StructVector *y );

/* struct_hig_innerprod.c */
HYPRE_Real hypre_StructHIGInnerProd ( hypre_StructVector *x , hypre_StructVector *y );

#ifdef __cplusplus
}
#endif

#endif