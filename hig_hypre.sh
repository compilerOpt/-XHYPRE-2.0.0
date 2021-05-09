#BHEADER**********************************************************************
#* 
#* XHYPRE 2020.10.29
#*
#* autor: Chuanying Li
#*
#EHEADER**********************************************************************

#!/bin/bash

dir=$(cd `dirname $0`; pwd)
hypre="/../hypre-2.11.2/src"
hypre_dir=$dir$hypre

#echo $dir
#echo $hypre_dir

cp -a $dir/hig_mv $hypre_dir

cp $dir/parcsr_mv/_hypre_parcsr_mv.h $hypre_dir/parcsr_mv/_hypre_parcsr_mv.h
cp $dir/parcsr_mv/Makefile $hypre_dir/parcsr_mv/Makefile
cp $dir/parcsr_mv/par_csr_matvec.c $hypre_dir/parcsr_mv/par_csr_matvec.c
cp $dir/parcsr_mv/par_vector.c $hypre_dir/parcsr_mv/par_vector.c

cp $dir/sstruct_mv/_hypre_sstruct_mv.h $hypre_dir/sstruct_mv/_hypre_sstruct_mv.h
cp $dir/sstruct_mv/Makefile $hypre_dir/sstruct_mv/Makefile
cp $dir/sstruct_mv/sstruct_innerprod.c $hypre_dir/sstruct_mv/sstruct_innerprod.c
cp $dir/sstruct_mv/sstruct_matvec.c $hypre_dir/sstruct_mv/sstruct_matvec.c

cp $dir/struct_ls/_hypre_struct_ls.h $hypre_dir/struct_ls/_hypre_struct_ls.h
cp $dir/struct_ls/Makefile $hypre_dir/struct_ls/Makefile
cp $dir/struct_ls/pcg_struct.c $hypre_dir/struct_ls/pcg_struct.c

cp $dir/src/Makefile $hypre_dir/Makefile

cp $dir/lib/Makefile $hypre_dir/lib/Makefile

cp $dir/FEI_mv--fei-hypre/Makefile $hypre_dir/FEI_mv/fei-hypre/Makefile

cp $dir/FEI_mv--femli/Makefile $hypre_dir/FEI_mv/femli/Makefile

cp $dir/IJ_mv/Makefile $hypre_dir/IJ_mv/Makefile

cp $dir/parcsr_block_mv/Makefile $hypre_dir/parcsr_block_mv/Makefile

cp $dir/parcsr_ls/Makefile $hypre_dir/parcsr_ls/Makefile

cp $dir/sstruct_ls/Makefile $hypre_dir/sstruct_ls/Makefile
