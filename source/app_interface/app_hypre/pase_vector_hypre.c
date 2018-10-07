#include "pase_config.h"
#include "pase_vector_hypre.h"
#include "_hypre_parcsr_mv.h"

#define DEBUG_PASE_VECTOR_HYPRE 1

//=============================================================================
//========================= HYPRE INTERFACE ===================================
//=============================================================================

/**
 * @brief 由向量创建向量
 */
void*
HYPRE_Vector_data_create_by_vector(void *x)
{
  HYPRE_ParVector  x_hypre      = (HYPRE_ParVector)x;
  MPI_Comm         comm         = hypre_ParVectorComm(x_hypre);
  PASE_INT         global_size  = hypre_ParVectorGlobalSize(x_hypre);
  PASE_INT        *partitioning = hypre_ParVectorPartitioning(x_hypre);

  HYPRE_ParVector  y = hypre_ParVectorCreate(comm, global_size, partitioning);
  HYPRE_ParVectorInitialize(y);
  hypre_ParVectorSetPartitioningOwner(y, 0);
  return (void*)y; 
}

/**
 * @brief 由矩阵创建向量
 */
void*
HYPRE_Vector_data_create_by_matrix(void *A)
{
  HYPRE_ParCSRMatrix A_hypre      = (HYPRE_ParCSRMatrix)A;
  MPI_Comm           comm         = hypre_ParCSRMatrixComm(A_hypre);
  PASE_INT           global_size  = hypre_ParCSRMatrixGlobalNumRows(A_hypre);
  PASE_INT          *partitioning = NULL;
#ifdef HYPRE_NO_GLOBAL_PARTITION
#if 1
  partitioning = hypre_CTAlloc(HYPRE_Int, 2);
  partitioning[0] = hypre_ParCSRMatrixFirstRowIndex(A_hypre);
  partitioning[1] = hypre_ParCSRMatrixLastRowIndex (A_hypre) + 1;
#else
  HYPRE_ParCSRMatrixGetRowPartitioning(A_hypre, &partitioning);
#endif
#else
  HYPRE_ParCSRMatrixGetRowPartitioning(A_hypre, &partitioning);
#endif

  HYPRE_ParVector y_hypre = hypre_ParVectorCreate(comm, global_size, partitioning);
  HYPRE_ParVectorInitialize(y_hypre);
  hypre_ParVectorSetPartitioningOwner(y_hypre, 1);
  return (void*)y_hypre;
}

/**
 * @brief 向量销毁
 */
void
HYPRE_Vector_data_destroy(void *x)
{
  HYPRE_ParVectorDestroy((HYPRE_ParVector)x);
  x = NULL;
}

/**
 * @brief 向量复制
 * @param src  输入参数, 被复制的向量
 * @param dst  输出参数, 复制得到的向量
 */
void
HYPRE_Vector_data_copy(void *src, void *dst)
{
  HYPRE_ParVectorCopy((HYPRE_ParVector)src, (HYPRE_ParVector)dst);
}

/**
 * @brief 设置向量元素为常数
 * @param x  输入/输出参数, 元素被设置为常数的向量
 * @param a  输入参数, 常数
 */
void
HYPRE_Vector_data_set_constant_value(void *x, PASE_SCALAR a)
{
  HYPRE_ParVectorSetConstantValues((HYPRE_ParVector)x, a);
}

/**
 * @brief 设置向量元素为随机数
 * @param x    输入/输出参数, 元素被设置为随机数的向量
 * @param seed 输入参数，随机数种子
 */
void
HYPRE_Vector_data_set_random_value(void *x, PASE_INT seed)
{
  HYPRE_ParVectorSetRandomValues((HYPRE_ParVector)x, seed);
}

/**
 * @brief 向量数乘 x = a * x
 */
void
HYPRE_Vector_data_scale(PASE_SCALAR a, void *x)
{
  HYPRE_ParVectorScale(a, (HYPRE_ParVector)x);
}

/**
 * @brief 向量校正 y = a * x + y
 */
void
HYPRE_Vector_data_axpy(PASE_SCALAR a, void *x, void *y)
{
  HYPRE_ParVectorAxpy(a, (HYPRE_ParVector)x, (HYPRE_ParVector)y);
}

/**
 * @brief 计算向量内积
 * @return 向量内积
 */
PASE_SCALAR
HYPRE_Vector_data_get_inner_product(void *x, void *y)
{
  PASE_SCALAR prod = 0.0;
  HYPRE_ParVectorInnerProd((HYPRE_ParVector)x, (HYPRE_ParVector)y, &prod);
  return prod;
}

/**
 * @brief 获得向量全局长度
 */
PASE_INT
HYPRE_Vector_data_get_global_nrow(void *x)
{
  return hypre_ParVectorGlobalSize((HYPRE_ParVector)x);
}

//=============================================================================
//=============================================================================
//=============================================================================

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_data_operator_create_by_hypre"
/**
 * @brief  基于 hypre 生成 PASE 向量操作集合
 * @return PASE_VECTOR_DATA_OPERATOR
 */
PASE_VECTOR_DATA_OPERATOR
PASE_Vector_data_operator_create_by_hypre()
{
    return PASE_Vector_data_operator_assign(HYPRE_Vector_data_create_by_vector,
                                            HYPRE_Vector_data_create_by_matrix,
                                            HYPRE_Vector_data_destroy,
                                            HYPRE_Vector_data_copy,
                                            HYPRE_Vector_data_set_constant_value,
                                            HYPRE_Vector_data_set_random_value,
                                            HYPRE_Vector_data_scale,
                                            HYPRE_Vector_data_axpy,
                                            HYPRE_Vector_data_get_inner_product,
                                            HYPRE_Vector_data_get_global_nrow);
} // end PASE_Vector_data_operator_create_by_hypre()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_create_by_hypre"
/**
 * @brief  基于 hypre 格式的向量数据生成 PASE 向量
 * @param  vector_data  输入参数, void * 指针, 指向向量数据
 * @return PASE_VECTOR
 */
PASE_VECTOR
PASE_Vector_create_by_hypre(void *vector_data)
{
  PASE_VECTOR_DATA_OPERATOR ops = PASE_Vector_data_operator_create_by_hypre();
  PASE_VECTOR x = PASE_Vector_assign(vector_data, ops);
  x->is_vector_data_owner = PASE_NO;
  x->data_interpreter     = DATA_INTERPRETER_HYPRE;
  PASE_Vector_data_operator_destroy(&ops);
  return x;
} // end PASE_Vector_create_by_hypre()