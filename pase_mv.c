/*
 * =====================================================================================
 *
 *       Filename:  pase_mv.c
 *
 *    Description:  PASE_ParCSRMatrix PASE_ParVetor各种操作
 *                  如果有任何内积操作的子函数, 都可以进行MPI的提升, 
 *                  即, 先做局部内积, 然后进行非阻塞归约求和, 之后进行其它操作
 *                  最后进行MPI_Wait, 并且应用得到的内积进行其它操作
 *
 *        Version:  1.0
 *        Created:  2017年09月08日 13时13分17秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LIYU
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <math.h>
#include "pase_mv.h"


/**
 * @brief 创建PASE并行矩阵
 *
 * @param comm 
 * @param N_H            aux_HH的阶数
 * @param block_size      aux_Hh这个MultiVector的向量个数
 * @param A_H            粗网格矩阵
 * @param A_h            细网格矩阵
 * @param P              P^T A_h P = A_H
 * @param u_h            要进行校正的向量组, 一共是block_size个
 * @param matrix         基于上述输入, 得到块矩阵
 *
 * @return               如果blockSize为0的情形也要考虑到
 *
 * 
 * 本函数主要是对PASE_ParCSRMatrix的成员赋值, 检查哪些是指针哪些需要malloc
 *
 */
PASE_Int PASE_ParCSRMatrixCreate( MPI_Comm comm , 
                                   PASE_Int block_size,
                                   HYPRE_ParCSRMatrix A_H, 
                                   HYPRE_ParCSRMatrix P,
                                   HYPRE_ParCSRMatrix A_h, 
				   HYPRE_ParVector*   u_h, 
				   PASE_ParCSRMatrix* matrix, 
				   HYPRE_ParVector    workspace_H, 
				   HYPRE_ParVector    workspace_h
				   )
{

   PASE_Int row, col;
   /* calloc并返回一个(pase_ParCSRMatrix *)的指针 */
   (*matrix) = hypre_CTAlloc(pase_ParCSRMatrix, 1);
   /* 是否需要为调用成员给出一个宏 */
   (*matrix)->comm = comm;
   (*matrix)->N_H = A_H->global_num_rows;
   (*matrix)->block_size = block_size;
   (*matrix)->A_H = A_H;
   (*matrix)->A_h = A_h;
   (*matrix)->P   = P;
   
   (*matrix)->aux_hH = hypre_CTAlloc(HYPRE_ParVector, block_size);
   (*matrix)->diag = 0;

   PASE_Int *partitioning;
   partitioning = hypre_ParVectorPartitioning(workspace_H);
   /* 创建并行向量 */
   for (col = 0; col < block_size; ++col)
   {
      (*matrix)->aux_hH[col] = hypre_ParVectorCreate(comm, (*matrix)->N_H, partitioning);
      hypre_ParVectorInitialize((*matrix)->aux_hH[col]);
      hypre_ParVectorOwnsPartitioning((*matrix)->aux_hH[col]) = 0;
   }

   (*matrix)->aux_Hh = (*matrix)->aux_hH;

   PASE_Int   num_nonzeros = block_size*block_size;
   (*matrix)->aux_hh = hypre_CSRMatrixCreate(block_size, block_size, num_nonzeros);
   hypre_CSRMatrixInitialize( (*matrix)->aux_hh );


   PASE_Int*  matrix_i = (*matrix)->aux_hh->i;
   /*第row行的非零元列号 matrix_j[ matrix_i[row] ]到matrix_j[ matrix_i[row+1] ] */
   PASE_Int*  matrix_j = (*matrix)->aux_hh->j;
   PASE_Real* matrix_data = (*matrix)->aux_hh->data;


   for (row = 0; row < block_size+1; ++row)
   {
      matrix_i[row] = row*block_size;
   }

   for ( row = 0; row < block_size; ++row)
   {
      for ( col = 0; col < block_size; ++col)
      {
	 matrix_j[col+row*block_size] = col;
      }
   }

//   for (row = 0; row < block_size; ++row)
//   {
//      /* y = alpha*A*x + beta*y */
//      hypre_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_h);
//      for ( col = 0; col < block_size; ++col)
//      {
//	 matrix_data[row*block_size+col] = hypre_ParVectorInnerProd(workspace_h, u_h[col]);
//      }
//
//      if (P==NULL)
//      {
//	 hypre_ParVectorCopy( workspace_h, (*matrix)->aux_hH[row] );
//      } 
//      else 
//      {
//	 hypre_ParCSRMatrixMatvecT(1.0, P, workspace_h, 0.0, (*matrix)->aux_hH[row]);
//      }
//   }

   if (u_h!=NULL)
   {
      /* MPI IMPROVEMENT */
      MPI_Status   status;
      MPI_Request *requests;
      requests = hypre_CTAlloc (MPI_Request, block_size);
      for (row = 0; row < block_size; ++row)
      {
	 hypre_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_h);
	 for ( col = row; col < block_size; ++col)
	 {
	    matrix_data[row*block_size+col] = hypre_SeqVectorInnerProd(
		  hypre_ParVectorLocalVector(workspace_h), hypre_ParVectorLocalVector(u_h[col]) );
	 }
	 MPI_Iallreduce( MPI_IN_PLACE, &matrix_data[row*block_size+row], block_size-row, MPI_DOUBLE, MPI_SUM, comm, &requests[row] );

	 if ( P == NULL )
	 {
	    hypre_ParVectorCopy(workspace_h, (*matrix)->aux_hH[row]);
	 } 
	 else 
	 {
	    hypre_ParCSRMatrixMatvecT(1.0, P, workspace_h, 0.0, (*matrix)->aux_hH[row]);
	 }
      }
      /* 从列开始 */
      for ( col = 0; col < block_size; ++col)
      {
	 MPI_Wait(&requests[col], &status);
	 for (row = col+1; row < block_size; ++row)
	 {
	    matrix_data[row*block_size+col] = matrix_data[col*block_size+row];
	 }
      }
      hypre_TFree(requests);
      /* MPI IMPROVEMENT */
   }
   return 0;
}


/**
 * @brief PASE矩阵基于block_size个u_h重新生成辅助部分
 * 
 * matrix->aux_hH[col] = P^T A_h u_h[row]
 * matrix->aux_hh  = u_h[col] A_h u_h[row]
 *
 * MPI IMPROVEMENT
 * 首先计算 workspace_h = A_h u_h[row] 
 * 对上三角部分的 matrix->aux_hh 进行局部内积计算 (由于是对称矩阵, 所以计算一半)
 * 计算完一行就进行全局非阻塞归约, (另外一种做法是全部算完之后再总的做一次归约)
 * 然后算matrix->aux_hH = P^T workspace_h
 *
 * 循环结束之后, 再对aux_hh的下三角部分进行赋值, 
 * 这里需要进行MPI_Wait, 然后再赋值
 *
 * @param comm
 * @param matrix
 * @param block_size
 * @param P
 * @param A_h
 * @param u_h
 * @param workspace_H
 * @param 
 *
 * @return 
 */
PASE_Int PASE_ParCSRMatrixSetAuxSpace( MPI_Comm      comm, 
				   PASE_ParCSRMatrix  matrix, 
                                   PASE_Int          block_size,
                                   HYPRE_ParCSRMatrix P,
                                   HYPRE_ParCSRMatrix A_h, 
				   HYPRE_ParVector*   u_h, 
				   PASE_Int          begin_idx,
				   HYPRE_ParVector    workspace_H, 
				   HYPRE_ParVector    workspace_h
				   )
{
/*  
   PASE_Int row, col;
   PASE_Real* matrix_data = matrix->aux_hh->data;
   for (row = 0; row < block_size; ++row)
   {
      hypre_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_h);
      for ( col = 0; col < block_size; ++col)
      {
	 matrix_data[row*block_size+col] = hypre_ParVectorInnerProd(workspace_h, u_h[col]);
      }
      if ( P == NULL )
      {
	 hypre_ParVectorCopy(workspace_h, matrix->aux_hH[row]);
      } 
      else 
      {
	 hypre_ParCSRMatrixMatvecT(1.0, P, workspace_h, 0.0, matrix->aux_hH[row]);
      }
   }
*/  
   /* MPI IMPROVEMENT */
   MPI_Status   status;
   MPI_Request *requests;
   requests = hypre_CTAlloc (MPI_Request, block_size);
   PASE_Int row, col, start;
   PASE_Real* matrix_data = matrix->aux_hh->data;
   for (row = 0; row < block_size; ++row)
   {
      /* TODO: 可以保存A_h u_h[row], row<begin_idx部分 */
      hypre_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_h);
      start = row>begin_idx?row:begin_idx;
      for ( col = start; col < block_size; ++col)
      {
	 matrix_data[row*block_size+col] = hypre_SeqVectorInnerProd(
	       hypre_ParVectorLocalVector(workspace_h), hypre_ParVectorLocalVector(u_h[col]) );
      }
      MPI_Iallreduce( MPI_IN_PLACE, &matrix_data[row*block_size+start], block_size-start, MPI_DOUBLE, MPI_SUM, comm, &requests[row] );

      if (row >= begin_idx)
      {
	 if ( P == NULL )
	 {
	    hypre_ParVectorCopy(workspace_h, matrix->aux_hH[row]);
	 } 
	 else 
	 {
	    hypre_ParCSRMatrixMatvecT(1.0, P, workspace_h, 0.0, matrix->aux_hH[row]);
	 }
      }
   }

   /* 从列开始 */
   for ( col = 0; col < block_size; ++col)
   {
      start = (col+1)>begin_idx?(col+1):begin_idx;
      MPI_Wait(&requests[col], &status);
      for (row = start; row < block_size; ++row)
      {
	 matrix_data[row*block_size+col] = matrix_data[col*block_size+row];
      }
   }
   hypre_TFree(requests);
   /* MPI IMPROVEMENT */

   return 0;
}


/**
 * @brief PASE_ParCSRMatrix销毁
 *
 * @param matrix
 *
 * @return 
 *
 * 将各个成员, 若是指针则置为NULL, 若有malloc则free
 *
 */
PASE_Int PASE_ParCSRMatrixDestroy( PASE_ParCSRMatrix matrix )
{

   int col;
   /* 是否需要为调用成员给出一个宏 */
   for (col = 0; col < matrix->block_size; ++col)
   {
      hypre_ParVectorDestroy(matrix->aux_hH[col]);
   }
   hypre_TFree(matrix->aux_hH);

   hypre_CSRMatrixDestroy(matrix->aux_hh);

   matrix->N_H = 0;
   matrix->block_size = 0;
   matrix->A_H = NULL;
   matrix->A_h = NULL;
   matrix->P   = NULL;
   matrix->aux_Hh = NULL;
   matrix->aux_hH = NULL;

   hypre_TFree(matrix);

   return 0;
}

PASE_Int PASE_ParCSRMatrixPrint( PASE_ParCSRMatrix matrix , const char *file_name )
{
   PASE_Int col, myid;
   MPI_Comm_rank(matrix->comm, &myid);
   char new_file_name[80];
   hypre_ParCSRMatrixPrint(matrix->A_H, file_name);
   hypre_sprintf(new_file_name,"%s.%s",file_name,"aux"); 
   if (myid==0)
   {
      hypre_CSRMatrixPrint(matrix->aux_hh, new_file_name);
   }
   for (col = 0; col < matrix->block_size; ++col)
   {
      hypre_sprintf(new_file_name,"%s.%d",new_file_name,col); 
      hypre_ParVectorPrint(matrix->aux_hH[col], new_file_name);
   }
   return 0;
}


PASE_Int PASE_ParVectorPrint( PASE_ParVector vector , const char *file_name )
{
   PASE_Int myid;
   MPI_Comm_rank(vector->comm, &myid);
   char new_file_name[80];
   hypre_ParVectorPrint(vector->b_H, file_name);
   hypre_sprintf(new_file_name,"%s.%s",file_name,"aux"); 
   if (myid==0)
   {
      hypre_SeqVectorPrint(vector->aux_h, new_file_name);
   }
   return 0;
}


PASE_Int PASE_ParVectorGetParVector( HYPRE_ParCSRMatrix P, PASE_Int block_size, 
      HYPRE_ParVector *vector_h, PASE_ParVector vector_Hh, HYPRE_ParVector vector )
{
   PASE_Int k;
   if ( P == NULL )
   {
      HYPRE_ParVectorCopy ( vector_Hh->b_H, vector );
      for (k = 0; k < block_size; ++k)
      {
	 /* y = a x + y */
	 HYPRE_ParVectorAxpy(vector_Hh->aux_h->data[k], vector_h[k], vector );
      }
   } 
   else 
   {
      HYPRE_ParCSRMatrixMatvec ( 1.0, P, vector_Hh->b_H, 0.0, vector );
      for (k = 0; k < block_size; ++k)
      {
	 /* y = a x + y */
	 HYPRE_ParVectorAxpy(vector_Hh->aux_h->data[k], vector_h[k], vector );
      }
   }
   return 0;
}


/**
 * @brief y = alpha A x + beta y
 *
 * @param alpha
 * @param A
 * @param x
 * @param beta
 * @param y
 *
 * @return 
 *
 * 这里x, y也是PASE_ParVector, 它们的成员b_H要与A的成员aux_HH的分布一致 
 * 判断矩阵和向量之间是否可以进行运算
 *
 * A_H    aux_hH    x->b_H      y->b_H          |   alpha A_H x->b_H    + alpha aux_hH x->aux_h  + beta y->b_H
 * aux_hH aux_hh    x->aux_h    y->aux_h        |   alpha aux_hH x->b_H + alpha aux_hh x->aux_h  + beta y->aux_h
 *
 * 为了提高并行效率, 先计算
 * A->aux_hH[col]与x->b_H进行内积: 
 * 每个进程进行局部向量的内积, 存储到tmp[col]里(PASE_Real数组)
 * 然后进行MPI_Iallreduce进行非阻塞全局归约
 * 之后, 进行与tmp无关的其它运算:
 * y->b_H = alpha A->A_H x->b_H + beta y->b_H
 * y->b_H = x->aux_h->data[col]) A->aux_hH[col] +  y->b_H
 *
 * y->aux_h = alpha A->aux_hh x->aux_h + beta y->aux_h
 * 
 * y->aux_h->data[col] += alpha * tmp;
 */
PASE_Int PASE_ParCSRMatrixMatvec ( PASE_Real alpha, PASE_ParCSRMatrix A, PASE_ParVector x, PASE_Real beta, PASE_ParVector y )
{

   PASE_Int col;
   PASE_Int block_size = A->block_size;

   if (A->diag == 0)
   {
      /* IMPROVEMENT for MPI */
      MPI_Comm    comm = hypre_ParVectorComm(x->b_H);
      MPI_Request request;
      MPI_Status  status;
      PASE_Real *inner_prod;
      inner_prod = hypre_CTAlloc (PASE_Real,  block_size);
      for (col = 0; col < block_size; ++col)
      {
	 inner_prod[col] = hypre_SeqVectorInnerProd(
	       hypre_ParVectorLocalVector(A->aux_hH[col]),  hypre_ParVectorLocalVector(x->b_H) );
      }
      MPI_Iallreduce(MPI_IN_PLACE, inner_prod, block_size, MPI_DOUBLE, MPI_SUM, comm, &request );
      /* IMPROVEMENT for MPI */

      /* y = alpha A x + beta y */
      hypre_ParCSRMatrixMatvec( alpha , A->A_H , x->b_H , beta , y->b_H );
      for (col = 0; col < block_size; ++col)
      {
	 /* y = alpha x + y */
	 hypre_ParVectorAxpy(alpha*(x->aux_h->data[col]), A->aux_hH[col], y->b_H);
      }
      hypre_CSRMatrixMatvec( alpha , A->aux_hh , x->aux_h , beta , y->aux_h );

      /* IMPROVEMENT for MPI */
      MPI_Wait(&request, &status);
      for (col = 0; col < block_size; ++col)
      {
	 //      y->aux_h->data[col] += alpha * hypre_ParVectorInnerProd(A->aux_hH[col], x->b_H);
	 y->aux_h->data[col] += alpha * inner_prod[col];
      }
      hypre_TFree(inner_prod);
      /* IMPROVEMENT for MPI */
   }
   else
   {
      hypre_ParCSRMatrixMatvec( alpha , A->A_H , x->b_H , beta , y->b_H );
      hypre_CSRMatrixMatvec( alpha , A->aux_hh , x->aux_h , beta , y->aux_h );
   }
   return 0;
}

/**
 * @brief y = alpha A^T x + beta y
 *
 * @param alpha
 * @param A
 * @param x
 * @param beta
 * @param y
 *
 * @return 
 *
 * 这里x也是PASE_ParVector
 * 判断矩阵和向量之间是否可以进行运算
 *
 */
PASE_Int PASE_ParCSRMatrixMatvecT( PASE_Real alpha, PASE_ParCSRMatrix A, PASE_ParVector x, PASE_Real beta, PASE_ParVector y )
{

   PASE_ParCSRMatrixMatvec ( alpha, A, x, beta, y );
   return 0;
}

PASE_Int PASE_ParVectorCopy(PASE_ParVector x, PASE_ParVector y)
{
   hypre_ParVectorCopy(x->b_H, y->b_H);
   hypre_SeqVectorCopy(x->aux_h, y->aux_h);
   return 0;
}

PASE_Int PASE_ParVectorInnerProd( PASE_ParVector x, PASE_ParVector y, PASE_Real *prod)
{
   /* TODO:可进行提升, 并行内积过程中计算aux_h的内积 */
   *prod = hypre_ParVectorInnerProd(x->b_H, y->b_H);
   *prod += hypre_SeqVectorInnerProd(x->aux_h, y->aux_h);
   return 0;
}

/* y = a x + y */
PASE_Int PASE_ParVectorAxpy( PASE_Real alpha , PASE_ParVector x , PASE_ParVector y )
{
   hypre_ParVectorAxpy(alpha, x->b_H, y->b_H); 
   hypre_SeqVectorAxpy(alpha, x->aux_h, y->aux_h); 
   return 0;
}


PASE_Int PASE_ParVectorSetConstantValues( PASE_ParVector v , PASE_Real value )
{
   hypre_ParVectorSetConstantValues( v->b_H , value );
   hypre_SeqVectorSetConstantValues( v->aux_h, value );
   return 0;
}

PASE_Int PASE_ParVectorSetRandomValues( PASE_ParVector v, PASE_Int seed )
{
   hypre_ParVectorSetRandomValues( v->b_H, seed );
   /* 这里因为每个进程是独立完成下述随机赋值, 但是由于seed重新给出, 
    * 所以每个进程得到的随机数是完全一样的, 
    * 这里seed+314的原因是与本进程上一条语句生成的随机数不同 */
   hypre_SeqVectorSetRandomValues( v->aux_h, seed+314 );
//   int i;
//   printf ( "SetRandomValues for PASE\n" );
//   for (i = 0; i<v->aux_h->size;++i)
//   {
//      printf ( "%f\n", v->aux_h->data[i] );
//   }
   return 0;
}





PASE_Int PASE_ParVectorScale ( PASE_Real alpha , PASE_ParVector y )
{
   hypre_ParVectorScale ( alpha , y->b_H );
   hypre_SeqVectorScale ( alpha , y->aux_h );
   return 0;
}

/**
 * @brief y = alpha A^T x + beta y 函数名的意思是强调向量的类型
 *
 * @param alpha
 * @param A           PASE_ParCSRMatrix
 * @param x           HYPRE_ParVector
 * @param beta
 * @param y           HYPRE_ParVector
 *
 * @return 
 *  
 * PASE并行块矩阵与HYPRE_Vector的矩阵向量运算, 会用在自己编写的PCG程序中, 以及调用Parpack时进行矩阵乘以向量
 * 如果可以利用HYPRE_ParVector进行PCG的编写, 则将HYPRE_Vector换做HYPRE_ParVector, 这样更对
 * 这里要判断矩阵向量是否可以进行运算.
 *
 * 大致算法:
 * 首先将x全收集(allgather)形成全局向量, 然后将之分为两个部分, 再进行运算, 多用指针
 * 一定要注意的是, 这里x, y的分布式存储和A的分布式存储完全不同.
 *
 */
PASE_Int PASE_ParCSRMatrixMatvec_HYPRE_ParVector ( PASE_Real alpha, PASE_ParCSRMatrix A, HYPRE_ParVector x, PASE_Real beta, HYPRE_ParVector y )
{

   return 0;
}

PASE_Int PASE_ParCSRMatrixMatvec_HYPRE_Vector ( PASE_Real alpha, PASE_ParCSRMatrix A, HYPRE_Vector x, PASE_Real beta, HYPRE_Vector y )
{

   return 0;
}


PASE_Int PASE_ParVectorCreate(    MPI_Comm comm , 
                                   PASE_Int N_H,
                                   PASE_Int block_size,
                                   HYPRE_ParVector    b_H, 
				   PASE_Int *partitioning, 
				   PASE_ParVector*    vector
				   )
{
   /* calloc并返回一个(pase_ParVector *)的指针 */
   (*vector) = hypre_CTAlloc(pase_ParVector, 1);
   /* 是否需要为调用成员给出一个宏 */
   (*vector)->comm = comm;
   if (b_H!=NULL)
   {
      (*vector)->N_H = b_H->global_size;
      (*vector)->b_H = b_H;
      (*vector)->owns_ParVector = 0;
   }
   else if (partitioning!=NULL) {
      (*vector)->N_H = N_H;
      (*vector)->b_H = hypre_ParVectorCreate( comm, N_H, partitioning );
      hypre_ParVectorInitialize((*vector)->b_H);
      hypre_ParVectorSetPartitioningOwner((*vector)->b_H,0);
      (*vector)->owns_ParVector = 1;
   }
   else {
      printf ( "Can't create PASE_ParVector without partitioning or HYPRE_ParVector.\n" );
      return -1;
   }

   (*vector)->block_size = block_size;
   (*vector)->aux_h = hypre_SeqVectorCreate(block_size);
   hypre_SeqVectorInitialize((*vector)->aux_h);

   return 0;
}

PASE_Int PASE_ParVectorDestroy( PASE_ParVector vector )
{
   hypre_SeqVectorDestroy(vector->aux_h);
   vector->aux_h = NULL;
   if (vector->owns_ParVector)
   {
      hypre_ParVectorDestroy(vector->b_H);
   }
   else {
      vector->b_H = NULL;
   }
   vector->block_size = 0;
   vector->N_H = 0;

   hypre_TFree(vector);
   return 0;
}




/**
 * @brief 
 *
 * TODO:是否可以对求内积部分进行MPI的改进
 *
 * 这里将pase_vertor的内积操作写成两部分, 其中对hypre部分按照
 * PASE_ParCSRMatrixSetAuxSpace中的改进一般, hypre_vector每个进程进行局部内积存入matrix_data
 * 同时要记录workspace_hH->aux_h与u_h[col]->aux_h的内积, 依然可以只计算上三角部分
 * 然后非阻塞进行全局归约求和, 之后进行与matrix_data无关的操作
 * 最后MPI_Wait, 之后加入workspace_hH->aux_h与u_h[col]->aux_h的内积, 并为下三角赋值
 *
 * @param comm
 * @param matrix
 * @param block_size
 * @param P
 * @param A_h
 * @param u_h
 * @param workspace_H
 * @param 
 *
 * @return 
 */
PASE_Int
PASE_ParCSRMatrixSetAuxSpaceByPASE_ParCSRMatrix( MPI_Comm comm , 
				   PASE_ParCSRMatrix   matrix, 
                                   PASE_Int           block_size,
                                   HYPRE_ParCSRMatrix  P,
                                   PASE_ParCSRMatrix   A_h, 
				   PASE_ParVector*     u_h, 
				   HYPRE_ParVector     workspace_H, 
				   PASE_ParVector      workspace_hH
				   )
{
   PASE_Int row, col;
   PASE_Real* matrix_data = matrix->aux_hh->data;
   for (row = 0; row < block_size; ++row)
   {
      /* y = alpha*A*x + beta*y */
      PASE_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_hH);
      for ( col = 0; col < block_size; ++col)
      {
	 PASE_ParVectorInnerProd(workspace_hH, u_h[col], &matrix_data[row*block_size+col] );
      }
      if ( P == NULL )
      {
	 hypre_ParVectorCopy(workspace_hH->b_H, matrix->aux_hH[row]);
      } 
      else
      {
	 hypre_ParCSRMatrixMatvecT(1.0, P, workspace_hH->b_H, 0.0, matrix->aux_hH[row]);
      }
   }
   return 0;
}


/**
 * @brief 以PASE_ParCSRMatrix为细空间矩阵进行粗空间PASE_ParCSRMatrix的生成
 *
 * A_h * u_h[i] 是PASE_ParVector, 它的b_H部分投影到粗空间生成A_H的aux_hH[i]部分
 * u_h[j] A_h * u_h[i] 生车A_H的aux_hh部分
 *
 * @param comm
 * @param block_size
 * @param P
 * @param A_h
 * @param u_h
 * @param matrix
 * @param workspace_H
 * @param workspace_hH
 *
 * @return 
 */
PASE_Int 
PASE_ParCSRMatrixCreateByPASE_ParCSRMatrix( MPI_Comm comm , 
                                   PASE_Int block_size,
                                   HYPRE_ParCSRMatrix A_H, 
                                   HYPRE_ParCSRMatrix P,
                                   PASE_ParCSRMatrix A_h, 
				   PASE_ParVector*   u_h, 
				   PASE_ParCSRMatrix* matrix, 
				   HYPRE_ParVector    workspace_H, 
				   PASE_ParVector    workspace_hH
				   )
{

   PASE_Int row, col;
   /* calloc并返回一个(pase_ParCSRMatrix *)的指针 */
   (*matrix) = hypre_CTAlloc(pase_ParCSRMatrix, 1);
   /* 是否需要为调用成员给出一个宏 */
   (*matrix)->comm = comm;
   (*matrix)->N_H = A_H->global_num_rows;
   (*matrix)->block_size = block_size;
   (*matrix)->A_H = A_H;
   /* 这个成员其实没有用 */
   (*matrix)->A_h = NULL;
   (*matrix)->P   = P;
   
   (*matrix)->aux_hH = hypre_CTAlloc(HYPRE_ParVector, block_size);


   PASE_Int *partitioning;
   partitioning = hypre_ParVectorPartitioning(workspace_H);
   /* 创建并行向量 */
   for (col = 0; col < block_size; ++col)
   {
      (*matrix)->aux_hH[col] = hypre_ParVectorCreate(comm, (*matrix)->N_H, partitioning);
      hypre_ParVectorInitialize((*matrix)->aux_hH[col]);
      hypre_ParVectorOwnsPartitioning((*matrix)->aux_hH[col]) = 0;
   }

   (*matrix)->aux_Hh = (*matrix)->aux_hH;

   PASE_Int   num_nonzeros = block_size*block_size;
   (*matrix)->aux_hh = hypre_CSRMatrixCreate(block_size, block_size, num_nonzeros);
   hypre_CSRMatrixInitialize( (*matrix)->aux_hh );


   PASE_Int*  matrix_i = (*matrix)->aux_hh->i;
   /*第row行的非零元列号 matrix_j[ matrix_i[row] ]到matrix_j[ matrix_i[row+1] ] */
   PASE_Int*  matrix_j = (*matrix)->aux_hh->j;
   PASE_Real* matrix_data = (*matrix)->aux_hh->data;


   for (row = 0; row < block_size+1; ++row)
   {
      matrix_i[row] = row*block_size;
   }

   for ( row = 0; row < block_size; ++row)
   {
      for ( col = 0; col < block_size; ++col)
      {
	 matrix_j[col+row*block_size] = col;
      }
   }

   /* TODO:同样可以进行MPI的提升 */
   for (row = 0; row < block_size; ++row)
   {
      /* y = alpha*A*x + beta*y */
      PASE_ParCSRMatrixMatvec(1.0, A_h, u_h[row], 0.0, workspace_hH);
//      PASE_ParCSRMatrixPrint(A_h, "workspace_A_h");
//      PASE_ParVectorPrint(u_h[row], "workspace_u_h");
//      PASE_ParVectorPrint(workspace_hH, "workspace_hH");
      for ( col = 0; col < block_size; ++col)
      {
	 PASE_ParVectorInnerProd(workspace_hH, u_h[col], &matrix_data[row*block_size+col] );
      }

      if (P==NULL)
      {
	 hypre_ParVectorCopy( workspace_hH->b_H, (*matrix)->aux_hH[row] );
      } 
      else 
      {
	 hypre_ParCSRMatrixMatvecT(1.0, P, workspace_hH->b_H, 0.0, (*matrix)->aux_hH[row]);
      }
   }

   return 0;
}


