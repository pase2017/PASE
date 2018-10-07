#include "pase_config.h"
#include "pase_matrix_hypre.h"
#include "_hypre_parcsr_mv.h"

#define DEBUG_PASE_MATRIX_HYPRE 1

//=============================================================================
//========================= HYPRE INTERFACE ===================================
//=============================================================================

/**
 * @brief 由矩阵创建矩阵
 */
void *
HYPRE_Matrix_data_create_by_matrix(void *A)
{
  HYPRE_ParCSRMatrix B = hypre_ParCSRMatrixCompleteClone((HYPRE_ParCSRMatrix)A);
  return (void*)B;
#if 0
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix)A, (void**) &parcsr_A);
  HYPRE_IJMatrix B;
  HYPRE_IJMatrixCreate(hypre_ParCSRMatrixComm         ((hypre_ParCSRMatrix*)parcsr_A),
                       hypre_ParCSRMatrixFirstRowIndex((hypre_ParCSRMatrix*)parcsr_A),
                       hypre_ParCSRMatrixLastRowIndex ((hypre_ParCSRMatrix*)parcsr_A),
                       hypre_ParCSRMatrixFirstRowIndex((hypre_ParCSRMatrix*)parcsr_A),
                       hypre_ParCSRMatrixLastRowIndex ((hypre_ParCSRMatrix*)parcsr_A),
                       &B);
  HYPRE_IJMatrixSetObjectType(B, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(B);
  return (void*)B;
#endif
}

/**
 * @brief 矩阵销毁
 */
void
HYPRE_Matrix_data_destroy(void *A)
{
  if(NULL != A) {
    HYPRE_ParCSRMatrixDestroy((HYPRE_ParCSRMatrix)A);
    A = NULL;
  }
#if 0
  if(NULL != A) {
    HYPRE_IJMatrixDestroy((HYPRE_IJMatrix)A);
    A = NULL;
  }
#endif
}

/**
 * @brief 矩阵复制
 * @param src  输入参数, 被复制的矩阵
 * @param dst  输出参数, 复制得到的矩阵
 */
void
HYPRE_Matrix_data_copy(void *src, void *dst)
{
  // 下面函数的第三个参数 1 代表 copy_data
  hypre_ParCSRMatrixCopy((HYPRE_ParCSRMatrix)src, (HYPRE_ParCSRMatrix)dst, 1);
#if 0
  HYPRE_ParCSRMatrix parcsr_src;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix)src, (void**) &parcsr_src);

  HYPRE_ParCSRMatrix parcsr_dst;
  HYPRE_IJMatrixGetObject((HYPRE_IJMatrix)dst, (void**) &parcsr_dst);
  // 下面函数的第三个参数 1 代表 copy_data
  hypre_ParCSRMatrixCopy(parcsr_src, parcsr_dst, 1);
#endif
}

/**
 * @brief 矩阵转置
 * @return 转置后的矩阵
 */
void *
HYPRE_Matrix_data_transpose(void *A)
{
  HYPRE_ParCSRMatrix AT;
  hypre_ParCSRMatrixTranspose(A, &AT, 1);
  return (void*)AT;
}

/**
 * @brief 矩阵矩阵相乘 C = A * B
 * @return 矩阵乘积 C
 */
void *
HYPRE_Matrix_data_multiply_matrix_matrix(void *A, void *B)
{
  HYPRE_ParCSRMatrix C = hypre_ParMatmul((HYPRE_ParCSRMatrix)A, (HYPRE_ParCSRMatrix)B);
  MPI_Comm comm = hypre_ParCSRMatrixComm((HYPRE_ParCSRMatrix)A);
  PASE_INT num_procs;
  MPI_Comm_size(comm, &num_procs);
  if(num_procs > 1) {
    hypre_MatvecCommPkgCreate(C);
  }
  return (void*)C;
}

/**
 * @brief 矩阵向量相乘 y = A * x
 */
void
HYPRE_Matrix_data_multiply_matrix_vector(void *A, void *x, void *y)
{
  hypre_ParCSRMatrixMatvec(1.0, (HYPRE_ParCSRMatrix)A, (HYPRE_ParVector)x, 0.0, (HYPRE_ParVector)y);
}

/**
 * @brief 矩阵向量相乘 y = a * A * x  + b * y
 */
void
HYPRE_Matrix_data_multiply_matrix_vector_general(PASE_SCALAR a, void *A, void *x, PASE_SCALAR b, void *y)
{
  hypre_ParCSRMatrixMatvec(a, (HYPRE_ParCSRMatrix)A, (HYPRE_ParVector)x, b, (HYPRE_ParVector)y);
}

/**
 * @brief 矩阵矩阵相乘 C = A^T * B
 * @return 矩阵乘积 C
 */
void *
HYPRE_Matrix_data_multiply_matrixT_matrix(void *A, void *B)
{
  HYPRE_ParCSRMatrix C = hypre_ParTMatmul((HYPRE_ParCSRMatrix)A, (HYPRE_ParCSRMatrix)B);
  MPI_Comm comm = hypre_ParCSRMatrixComm((HYPRE_ParCSRMatrix)A);
  PASE_INT num_procs;
  MPI_Comm_size(comm, &num_procs);
  if(num_procs > 1) {
    hypre_MatvecCommPkgCreate(C);
  }
  return (void*)C;
}

/**
 * @brief 矩阵向量相乘 y = A^T * x
 */
void
HYPRE_Matrix_data_multiply_matrixT_vector(void *A, void *x, void *y)
{
  hypre_ParCSRMatrixMatvecT(1.0, (HYPRE_ParCSRMatrix)A, (HYPRE_ParVector)x, 0.0, (HYPRE_ParVector)y);
}

/**
 * @brief 矩阵向量相乘 y = a * A^T * x  + b * y
 */
void
HYPRE_Matrix_data_multiply_matrixT_vector_general(PASE_SCALAR a, void *A, void *x, PASE_SCALAR b, void *y)
{
  hypre_ParCSRMatrixMatvecT(a, (HYPRE_ParCSRMatrix)A, (HYPRE_ParVector)x, b, (HYPRE_ParVector)y);
}

/**
 * @brief 获得矩阵全局行数
 */
PASE_INT
HYPRE_Matrix_data_get_global_nrow(void *A)
{
  return hypre_ParCSRMatrixGlobalNumRows((HYPRE_ParCSRMatrix)A);
}

/**
 * @brief 获得矩阵全局列数
 */
PASE_INT
HYPRE_Matrix_data_get_global_ncol(void *A)
{
  return hypre_ParCSRMatrixGlobalNumCols((HYPRE_ParCSRMatrix)A);
}

/**
 * @brief 获得矩阵通信器
 */
MPI_Comm
HYPRE_Matrix_data_get_mpi_comm(void *A)
{
  return hypre_ParCSRMatrixComm((HYPRE_ParCSRMatrix)A);
}

//=============================================================================
//=============================================================================
//=============================================================================

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_data_operator_create_by_hypre"
/**
 * @brief 基于 hypre 生成 PASE 矩阵操作集合
 */
PASE_MATRIX_DATA_OPERATOR
PASE_Matrix_data_operator_create_by_hypre()
{
  return PASE_Matrix_data_operator_assign(HYPRE_Matrix_data_create_by_matrix,
                                          HYPRE_Matrix_data_destroy,
                                          HYPRE_Matrix_data_copy,
                                          HYPRE_Matrix_data_transpose,
                                          HYPRE_Matrix_data_multiply_matrix_matrix,
                                          HYPRE_Matrix_data_multiply_matrix_vector,
                                          HYPRE_Matrix_data_multiply_matrix_vector_general,
                                          HYPRE_Matrix_data_multiply_matrixT_matrix,
                                          HYPRE_Matrix_data_multiply_matrixT_vector,
                                          HYPRE_Matrix_data_multiply_matrixT_vector_general,
                                          HYPRE_Matrix_data_get_global_nrow,
                                          HYPRE_Matrix_data_get_global_ncol,
                                          HYPRE_Matrix_data_get_mpi_comm);
} // PASE_Matrix_data_operator_create_by_hypre()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_create_by_hypre"
/**
 * @brief 基于 hypre 格式的矩阵数据生成 PASE 矩阵
 */
PASE_MATRIX
PASE_Matrix_create_by_hypre(void *matrix_data)
{
  PASE_MATRIX_DATA_OPERATOR ops = PASE_Matrix_data_operator_create_by_hypre();
  PASE_MATRIX A = PASE_Matrix_assign(matrix_data, ops);
  A->is_matrix_data_owner = PASE_NO;
  A->data_interpreter     = DATA_INTERPRETER_HYPRE;
  PASE_Matrix_data_operator_destroy(&ops);
  return A;
} // end PASE_Matrix_create_by_hypre()
