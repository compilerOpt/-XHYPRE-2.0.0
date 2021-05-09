/*BHEADER**********************************************************************
 * 
 * XHYPRE 2020.10.29
 *
 * autor: Chuanying Li
 *
 ***********************************************************************EHEADER*/


此tar包依赖于hypre库。

hypre库版本：hypre-2.11.2

作用为：将原始hypre库转换为高精度hypre库。

使用步骤：
首先将此tar包与hypre-2.11.2.tar.gz上传至同一文件夹下；
然后分别解压缩；
进入解压缩后的hig_hypre文件夹，运行脚本hig_hypre.sh，运行语句为：
                 chmod +x hig_hypre.sh
                 ./hig_hypre.sh
此时原始hypre就变成了高精度的hypre库。

按照原始hypre库编译运行的方式使用即可。