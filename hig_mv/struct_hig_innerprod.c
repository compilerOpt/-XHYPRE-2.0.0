/*BHEADER**********************************************************************
 * 
 * XHYPRE 2020.10.29
 *
 * autor: Chuanying Li
 *
 ***********************************************************************EHEADER*/

/******************************************************************************
 *
 * Structured inner product routine
 *
 *****************************************************************************/
 
#include "hig_mv.h"
#include "_hypre_struct_mv.h"

#ifndef SPLITTER
#define SPLITTER 134217729.0               // = 2^27 + 1
#endif

/*--------------------------------------------------------------------------
 * hypre_HIGStructInnerProd
 *--------------------------------------------------------------------------*/
 

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
	
HYPRE_Real
hypre_StructHIGInnerProd( hypre_StructVector *x,
                       hypre_StructVector *y )
{
   HYPRE_Real       final_innerprod_result;
   HYPRE_Real       local_result;
   HYPRE_Real       process_result;
                   
   hypre_Box       *x_data_box;
   hypre_Box       *y_data_box;
                   
   HYPRE_Int        xi;
   HYPRE_Int        yi;
                   
   HYPRE_Complex   *xp;
   HYPRE_Complex   *yp;
                   
   hypre_BoxArray  *boxes;
   hypre_Box       *box;
   hypre_Index      loop_size;
   hypre_IndexRef   start;
   hypre_Index      unit_stride;
                   
   HYPRE_Int        i,t;
   
   HYPRE_Complex    tempx;  

   local_result = 0.0;
   process_result = 0.0;
   
   dd_real temp1;   
   tempx = 0.0;
   dd_real temps= {0.0,0.0};

   hypre_SetIndex(unit_stride, 1);

   boxes = hypre_StructGridBoxes(hypre_StructVectorGrid(y));
   hypre_ForBoxI(i, boxes)
   {
      box   = hypre_BoxArrayBox(boxes, i);
      start = hypre_BoxIMin(box);

      x_data_box = hypre_BoxArrayBox(hypre_StructVectorDataSpace(x), i);
      y_data_box = hypre_BoxArrayBox(hypre_StructVectorDataSpace(y), i);

      xp = hypre_StructVectorBoxData(x, i);
      yp = hypre_StructVectorBoxData(y, i);

      hypre_BoxGetSize(box, loop_size);

      hypre_BoxLoop2Begin(hypre_StructVectorNDim(x), loop_size,
                          x_data_box, start, unit_stride, xi,
                          y_data_box, start, unit_stride, yi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,xi,yi) reduction(+:local_result) HYPRE_SMP_SCHEDULE
#endif
      
      hypre_BoxLoop2For(xi, yi)
      {
		  temp1 = two_prod(xp[xi], hypre_conj(yp[yi]));
		  temps = two_sum(temps.H, temp1.H);
		  tempx = tempx + (temp1.L + temps.L) ;
      }
      hypre_BoxLoop2End(xi, yi);
	  tempx = tempx + temps.H;
	  local_result += tempx;
   }
   process_result = local_result;

   hypre_MPI_Allreduce(&process_result, &final_innerprod_result, 1,
                       HYPRE_MPI_REAL, hypre_MPI_SUM, hypre_StructVectorComm(x));

   hypre_IncFLOPCount(2*hypre_StructVectorGlobalSize(x));

   return final_innerprod_result;
}

