/*BHEADER**********************************************************************
 * 
 * XHYPRE 2020.10.29
 *
 * autor: Chuanying Li
 *
 ***********************************************************************EHEADER*/
 

/******************************************************************************
 *
 * Member functions for hypre_Vector class.
 *
 *****************************************************************************/

#include "hig_mv.h"
#include "seq_mv.h"
#include <assert.h>

#ifndef SPLITTER
#define SPLITTER 134217729.0               // = 2^27 + 1
#endif

typedef struct doubledouble{
    double H;
    double L;
}dd_real;
	
/* Computes fl(a+b) and err(a+b).  */
static dd_real two_sum(double a, double b) {
    dd_real res;
	double bb;
    res.H = a + b;
    bb = res.H - a;
    res.L = (a - (res.H - bb)) + (b - bb);
    return res;
}
	
/* Computes fl(a*b) and err(a*b) by using FMA. */
static dd_real two_prod(double a, double b) {
    dd_real res;
    res.H = a * b;
    res.L = fma(a,b,-res.H);
    return res;
}


/*--------------------------------------------------------------------------
 * hypre_HIGSeqVectorInnerProd
 *--------------------------------------------------------------------------*/

HYPRE_Real   hypre_HIGSeqVectorInnerProd( hypre_Vector *x,
                                       hypre_Vector *y )
{
#ifdef HYPRE_PROFILE
   hypre_profile_times[HYPRE_TIMER_ID_BLAS1] -= hypre_MPI_Wtime();
#endif

   HYPRE_Complex *x_data = hypre_VectorData(x);
   HYPRE_Complex *y_data = hypre_VectorData(y);
   HYPRE_Int      size   = hypre_VectorSize(x);
           
   HYPRE_Int      i;

   HYPRE_Real     result = 0.0;
   
   dd_real temp1;

   size *=hypre_VectorNumVectors(x);
   

#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(i) reduction(+:result) HYPRE_SMP_SCHEDULE
#endif
   dd_real temps= {0.0,0.0};
   for (i = 0; i < size; i++)
   {
	   temp1 = two_prod(hypre_conj(y_data[i]), x_data[i]);
	   temps = two_sum(temps.H, temp1.H);
	   result = result + (temp1.L + temps.L) ;
   }
   result = result + temps.H;

#ifdef HYPRE_PROFILE
   hypre_profile_times[HYPRE_TIMER_ID_BLAS1] += hypre_MPI_Wtime();
#endif
   return result;
}
/* XHYPRE 2020.11.29 */

