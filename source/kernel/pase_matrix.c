#include "pase_config.h"
#include "pase_matrix.h"

#define DEBUG_PASE_MATRIX 1

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_data_operator_assign"
/**
 * @brief 指定矩阵操作集合
 */
PASE_MATRIX_DATA_OPERATOR 
PASE_Matrix_data_operator_assign
    (void * (*create_by_matrix)               (void *A),
     void   (*destroy)                        (void *A),
     void   (*copy)                           (void *A, void *B),
     void * (*transpose)                      (void *A),
     void * (*multiply_matrix_matrix)         (void *A, void *B),
     void   (*multiply_matrix_vector)         (void *A, void *x, void *y),
     void   (*multiply_matrix_vector_general) (PASE_SCALAR a, void *A, void *x, PASE_SCALAR b, void *y),
     void * (*multiply_matrixT_matrix)        (void *A, void *B),
     void   (*multiply_matrixT_vector)        (void *A, void *x, void *y),
     void   (*multiply_matrixT_vector_general)(PASE_SCALAR a, void *A, void *x, PASE_SCALAR b, void *y),
     PASE_INT (*get_global_nrow)              (void *A),
     PASE_INT (*get_global_ncol)              (void *A),
     MPI_Comm (*get_mpi_comm)                 (void *A))
{
  PASE_MATRIX_DATA_OPERATOR ops 
      = (PASE_MATRIX_DATA_OPERATOR)PASE_Malloc(sizeof(PASE_MATRIX_DATA_OPERATOR_PRIVATE));

  ops->create_by_matrix                = create_by_matrix;
  ops->destroy                         = destroy;
  ops->copy                            = copy;
  ops->transpose                       = transpose;
  ops->multiply_matrix_matrix          = multiply_matrix_matrix;
  ops->multiply_matrix_vector          = multiply_matrix_vector;
  ops->multiply_matrix_vector_general  = multiply_matrix_vector_general;
  ops->multiply_matrixT_matrix         = multiply_matrixT_matrix;
  ops->multiply_matrixT_vector         = multiply_matrixT_vector;
  ops->multiply_matrixT_vector_general = multiply_matrixT_vector_general;
  ops->get_global_nrow                 = get_global_nrow;
  ops->get_global_ncol                 = get_global_ncol;
  ops->get_mpi_comm                    = get_mpi_comm;

  return ops;
} // PASE_Matrix_data_operator_assign()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_data_operator_destroy"
/**
 * @brief 销毁矩阵操作集合
 */
void
PASE_Matrix_data_operator_destroy(PASE_MATRIX_DATA_OPERATOR *ops)
{
  PASE_Free(*ops);
} // end PASE_Matrix_data_operator_destroy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_assign"
/**
 * @brief 将矩阵数据与操作集合关联至矩阵
 */
PASE_MATRIX 
PASE_Matrix_assign(void *matrix_data, PASE_MATRIX_DATA_OPERATOR ops)
{
#if DEBUG_PASE_MATRIX
  if(NULL == matrix_data) {
    PASE_Error(__FUNCT__": Cannot create PASE MATRIX without matrix data.\n");
  }
  if(NULL == ops) {
    PASE_Error(__FUNCT__": Cannot create PASE MATRIX without matrix data operator.\n");
  }
  if(NULL == ops->get_global_nrow) {
    PASE_Error(__FUNCT__": ops->get_global_nrow is not assigned.\n");
  }
  if(NULL == ops->get_global_ncol) {
    PASE_Error(__FUNCT__": ops->get_global_ncol is not assigned.\n");
  }
#endif
  
  PASE_MATRIX A = (PASE_MATRIX)PASE_Malloc(sizeof(PASE_MATRIX_PRIVATE));
  A->ops = (PASE_MATRIX_DATA_OPERATOR)PASE_Malloc(sizeof(PASE_MATRIX_DATA_OPERATOR_PRIVATE));
  *(A->ops)               = *ops;
  A->matrix_data          = matrix_data;
  A->is_matrix_data_owner = PASE_NO;
  A->data_interpreter     = DATA_INTERPRETER_USER;
  A->global_nrow = A->ops->get_global_nrow(A->matrix_data);
  A->global_ncol = A->ops->get_global_ncol(A->matrix_data);
  return A;
} // end PASE_Matrix_assign()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_destroy"
/**
 * @brief 销毁矩阵
 */
void 
PASE_Matrix_destroy(PASE_MATRIX *A)
{
  if(NULL == *A) return;

#if DEBUG_PASE_MATRIX
  if((PASE_YES != (*A)->is_matrix_data_owner) &&
     (PASE_NO  != (*A)->is_matrix_data_owner)) {
    PASE_Error(__FUNCT__": Cannot decide whether the owner of matrix is.");
  }
  if(PASE_YES == (*A)->is_matrix_data_owner) {
    if(NULL == (*A)->ops->destroy) {
      PASE_Error(__FUNCT__": Operator ops->destroy is not assigned.\n");
    }
  }
#endif

  if(PASE_YES == (*A)->is_matrix_data_owner) {
    (*A)->ops->destroy((*A)->matrix_data);
  }
  PASE_Matrix_data_operator_destroy(&((*A)->ops));
  PASE_Free(*A);
} // end PASE_Matrix_destroy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_copy"
/**
 * @brief 矩阵复制
 * @note 1. 该函数不负责分配内存, 因此 src 和 dst 的 matrix_data 不能为空
 *       2. 仅用于 src 和 dst 具有相同稀疏结构的情况
 */
void 
PASE_Matrix_copy(PASE_MATRIX src, PASE_MATRIX dst)
{
#if DEBUG_PASE_MATRIX
  if((NULL == src) || (NULL == dst)) {
    PASE_Error(__FUNCT__": Neither the two matrices can be empty.\n");
  }
  if((NULL == src->matrix_data) || (NULL == dst->matrix_data)) {
    PASE_Error(__FUNCT__": Cannot copy matrix if matrix data of either src or dst is NULL.\n");
  }
  if(dst->data_interpreter != src->data_interpreter) {
    PASE_Error(__FUNCT__": Cannot copy a matrix to a different interpreter.\n");
  }
  if(NULL == src->ops->copy) {
    PASE_Error(__FUNCT__": Operator ops->copy is not assigned.\n");
  }
  if((src->global_nrow != dst->global_nrow) || (src->global_ncol != dst->global_ncol)) {
    PASE_Error(__FUNCT__": Cannot copy matrix if dimensions are not matched.\n");
  }
#endif

  src->ops->copy(src->matrix_data, dst->matrix_data);
} // end PASE_Matrix_copy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_get_global_nrow"
/**
 * @brief 获取矩阵全局行数
 */
PASE_INT
PASE_Matrix_get_global_nrow(PASE_MATRIX A)
{
#if DEBUG_PASE_MATRIX
  if(NULL == A) {
    PASE_Error(__FUNCT__": Matrix cannot be empty.\n");
  }
  if(NULL == A->ops->get_global_nrow) {
    PASE_Error(__FUNCT__": Operator ops->get_global_nrow is not assigned.\n");
  }
#endif
  return A->ops->get_global_nrow(A->matrix_data);
} // end PASE_Matrix_get_global_nrow()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_get_global_ncol"
/**
 * @brief 获取矩阵全局行数
 */
PASE_INT
PASE_Matrix_get_global_ncol(PASE_MATRIX A)
{
#if DEBUG_PASE_MATRIX
  if(NULL == A) {
    PASE_Error(__FUNCT__": Matrix cannot be empty.\n");
  }
  if(NULL == A->ops->get_global_ncol) {
    PASE_Error(__FUNCT__": Operator ops->get_global_ncol is not assigned.\n");
  }
#endif
  return A->ops->get_global_ncol(A->matrix_data);
} // end PASE_Matrix_get_global_ncol()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_get_mpi_comm"
/**
 * @brief 获取矩阵通信器
 */
MPI_Comm 
PASE_Matrix_get_mpi_comm(PASE_MATRIX A)
{
#if DEBUG_PASE_MATRIX
  if(NULL == A) {
    PASE_Error(__FUNCT__": Matrix cannot be empty.\n");
  }
  if(NULL == A->ops->get_mpi_comm) {
    PASE_Error(__FUNCT__": Operator ops->get_mpi_comm is not assigned.\n");
  }
#endif
  return A->ops->get_mpi_comm(A->matrix_data);
} // end PASE_Matrix_get_mpi_comm()
