/*BHEADER**********************************************************************
 * 
 * XHYPRE 2020.10.29
 *
 * autor: Chuanying Li
 *
 ***********************************************************************EHEADER*/

/******************************************************************************
 *
 * Structured matrix-vector multiply routine
 *
 *****************************************************************************/

#include "hig_mv.h"
#include "_hypre_struct_mv.h"

/* this currently cannot be greater than 7 */
#ifdef MAX_DEPTH
#undef MAX_DEPTH
#endif
#define MAX_DEPTH 7

#ifndef SPLITTER
#define SPLITTER 134217729.0               // = 2^27 + 1
#endif

/*--------------------------------------------------------------------------
 * hypre_StructMatvecData data structure
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

typedef struct
{
   hypre_StructMatrix  *A;
   hypre_StructVector  *x;
   hypre_ComputePkg    *compute_pkg;

} hypre_StructMatvecData;



/*--------------------------------------------------------------------------
 * hypre_StructHIGMatvecCompute
 *--------------------------------------------------------------------------*/


HYPRE_Int
hypre_StructHIGMatvecCompute( void               *matvec_vdata,
                           HYPRE_Complex       alpha,
                           hypre_StructMatrix *A,
                           hypre_StructVector *x,
                           HYPRE_Complex       beta,
                           hypre_StructVector *y            )
{
	hypre_StructMatvecData  *matvec_data = (hypre_StructMatvecData  *)matvec_vdata;
                          
   hypre_ComputePkg        *compute_pkg;
                          
   hypre_CommHandle        *comm_handle;
                          
   hypre_BoxArrayArray     *compute_box_aa;
   hypre_Box               *y_data_box;
                          
   HYPRE_Int                yi;
                          
   HYPRE_Complex           *xp;
   HYPRE_Complex           *yp;
                          
   hypre_BoxArray          *boxes;
   hypre_Box               *box;
   hypre_Index              loop_size;
   hypre_IndexRef           start;
   hypre_IndexRef           stride;
                          
   HYPRE_Int                constant_coefficient;

   HYPRE_Complex            temp;
   HYPRE_Int                compute_i, i;

   hypre_StructVector      *x_tmp = NULL;
   
   HYPRE_Complex    tempx;  
   dd_real temp1;  
   tempx = 0.0;
   dd_real temps= {0.0,0.0};
   

   /*-----------------------------------------------------------------------
    * Initialize some things
    *-----------------------------------------------------------------------*/

   constant_coefficient = hypre_StructMatrixConstantCoefficient(A);
   if (constant_coefficient) hypre_StructVectorClearBoundGhostValues(x, 0);

   compute_pkg = (matvec_data -> compute_pkg);

   stride = hypre_ComputePkgStride(compute_pkg);

   /*-----------------------------------------------------------------------
    * Do (alpha == 0.0) computation
    *-----------------------------------------------------------------------*/

   if (alpha == 0.0)
   {
      boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(A));
      hypre_ForBoxI(i, boxes)
      {
         box   = hypre_BoxArrayBox(boxes, i);
         start = hypre_BoxIMin(box);

         y_data_box = hypre_BoxArrayBox(hypre_StructVectorDataSpace(y), i);
         yp = hypre_StructVectorBoxData(y, i);

         hypre_BoxGetSize(box, loop_size);

         hypre_BoxLoop1Begin(hypre_StructVectorNDim(x), loop_size,
                             y_data_box, start, stride, yi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,yi) HYPRE_SMP_SCHEDULE
#endif
         hypre_BoxLoop1For(yi)
         {
            yp[yi] *= beta;
         }
         hypre_BoxLoop1End(yi);
      }

      return hypre_error_flag;
   }

   if (x == y)
   {
      x_tmp = hypre_StructVectorClone(y);
      x = x_tmp;
   }
   /*-----------------------------------------------------------------------
    * Do (alpha != 0.0) computation
    *-----------------------------------------------------------------------*/

   for (compute_i = 0; compute_i < 2; compute_i++)
   {
      switch(compute_i)
      {
         case 0:
         {
            xp = hypre_StructVectorData(x);
            hypre_InitializeIndtComputations(compute_pkg, xp, &comm_handle);
            compute_box_aa = hypre_ComputePkgIndtBoxes(compute_pkg);

            /*--------------------------------------------------------------
             * initialize y= (beta/alpha)*y normally (where everything
             * is multiplied by alpha at the end),
             * beta*y for constant coefficient (where only Ax gets multiplied by alpha)
             *--------------------------------------------------------------*/

            if ( constant_coefficient==1 )
            {
               temp = beta;
            }
            else
            {
               temp = beta / alpha;
            }
            if (temp != 1.0)
            {
               boxes = hypre_StructGridBoxes(hypre_StructMatrixGrid(A));
               hypre_ForBoxI(i, boxes)
               {
                  box   = hypre_BoxArrayBox(boxes, i);
                  start = hypre_BoxIMin(box);

                  y_data_box =
                     hypre_BoxArrayBox(hypre_StructVectorDataSpace(y), i);
                  yp = hypre_StructVectorBoxData(y, i);

                  if (temp == 0.0)
                  {
                     hypre_BoxGetSize(box, loop_size);

                     hypre_BoxLoop1Begin(hypre_StructVectorNDim(x), loop_size,
                                         y_data_box, start, stride, yi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,yi) HYPRE_SMP_SCHEDULE
#endif
                     hypre_BoxLoop1For(yi)
                     {
                        yp[yi] = 0.0;
                     }
                     hypre_BoxLoop1End(yi);
                  }
                  else
                  {
                     hypre_BoxGetSize(box, loop_size);

                     hypre_BoxLoop1Begin(hypre_StructVectorNDim(x), loop_size,
                                         y_data_box, start, stride, yi);
#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(HYPRE_BOX_PRIVATE,yi) HYPRE_SMP_SCHEDULE
#endif
                     
                     hypre_BoxLoop1For(yi)
                     {
						temp1 = two_prod(yp[yi], temp);
		                yp[yi] = temp1.L + temp1.H ; 
                     }
                     hypre_BoxLoop1End(yi);
                  }
               }
            }
         }
         break;

         case 1:
         {
            hypre_FinalizeIndtComputations(comm_handle);
            compute_box_aa = hypre_ComputePkgDeptBoxes(compute_pkg);
         }
         break;
      }

      /*--------------------------------------------------------------------
       * y += A*x
       *--------------------------------------------------------------------*/

      switch( constant_coefficient )
      {
         case 0:
         {
            hypre_StructMatvecCC0( alpha, A, x, y, compute_box_aa, stride );
            break;
         }
         case 1:
         {
            hypre_StructMatvecCC1( alpha, A, x, y, compute_box_aa, stride );
            break;
         }
         case 2:
         {
            hypre_StructMatvecCC2( alpha, A, x, y, compute_box_aa, stride );
            break;
         }
      }

   }

   if (x_tmp)
   {
      hypre_StructVectorDestroy(x_tmp);
      x = y;
   }
   
   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_StructHIGMatvec
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_StructHIGMatvec( HYPRE_Complex       alpha,
                    hypre_StructMatrix *A,
                    hypre_StructVector *x,
                    HYPRE_Complex       beta,
                    hypre_StructVector *y     )
{
   void *matvec_data;

   matvec_data = hypre_StructMatvecCreate();
   hypre_StructMatvecSetup(matvec_data, A, x);
   hypre_StructHIGMatvecCompute(matvec_data, alpha, A, x, beta, y);
   hypre_StructMatvecDestroy(matvec_data);

   return hypre_error_flag;
}
