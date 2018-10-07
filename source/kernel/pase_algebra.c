#include <math.h>
#include <assert.h>
#include "pase_algebra.h"

#define DEBUG_PASE_ALGEBRA 1

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_transpose"
/**
 * @brief 矩阵转置
 */
PASE_MATRIX 
PASE_Matrix_transpose(PASE_MATRIX A)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == A) {
    PASE_Error(__FUNCT__": Matrix cannot be empty.\n");
  }
  if(NULL == A->ops) {
    PASE_Error(__FUNCT__": A->ops is not assigned.\n");
  }
  if(NULL == A->ops->transpose) {
    PASE_Error(__FUNCT__": A->ops->transpose is not assigned.\n");
  }
#endif

  void *matrix_data        = A->ops->transpose(A->matrix_data);
  PASE_MATRIX AT           = PASE_Matrix_assign(matrix_data, A->ops); 
  AT->is_matrix_data_owner = PASE_YES;
  AT->data_interpreter     = A->data_interpreter;
  return AT;
} // end PASE_Matrix_transpose()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_multiply_matrix"
/**
 * @brief 矩阵与矩阵相乘 A * B
 */
PASE_MATRIX 
PASE_Matrix_multiply_matrix(PASE_MATRIX A, PASE_MATRIX B)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == B)) {
    PASE_Error(__FUNCT__": Neither the two matrices can be empty.\n");
  }
  if(A->global_ncol != B->global_nrow) {
    PASE_Error(__FUNCT__": Matrix dimensions are not match.\n");
  }
  if(A->data_interpreter != B->data_interpreter) {
    PASE_Error(__FUNCT__": Matrix data_interpreters are not match.\n");
  }
  if(NULL == A->ops->multiply_matrix_matrix) {
    PASE_Error(__FUNCT__": Operator ops->multiply_matrix_matrix is not assigned.\n");
  }
#endif

  void *matrix_data = A->ops->multiply_matrix_matrix(A->matrix_data, B->matrix_data); 
  PASE_MATRIX C = PASE_Matrix_assign(matrix_data, A->ops);
  C->is_matrix_data_owner = PASE_YES;
  C->data_interpreter     = A->data_interpreter;
  return C;
} // end PASE_Matrix_multiply_matrix()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_multiply_vector"
/**
 * @brief 矩阵向量相乘 y = A * x (其中 y 不能与 x 相同)
 */
void 
PASE_Matrix_multiply_vector(PASE_MATRIX A, PASE_VECTOR x, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Matrix and vectors cannot be empty.\n");
  }
  if((A->global_ncol != x->global_nrow) ||
     (A->global_nrow != y->global_nrow)) {
    PASE_Error(__FUNCT__": The dimensions of matrix and vector are not matched.\n");
  }
  if((NULL == A->ops->multiply_matrix_vector) &&
     (NULL == A->ops->multiply_matrix_vector_general)) {
    PASE_Error(__FUNCT__": Neither ops->multiply_matrix_vector nor "
                        "ops->multiply_matrix_vector_general is assigned.\n");
  }
  if(!(A->data_interpreter == x->data_interpreter && A->data_interpreter == y->data_interpreter)) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if(y == x) {
    PASE_Error(__FUNCT__": Vectors x and y cannot be the same.\n");
  }
#endif
  
  if(NULL != A->ops->multiply_matrix_vector) {
    A->ops->multiply_matrix_vector(A->matrix_data, x->vector_data, y->vector_data);
  } else {
    PASE_Matrix_multiply_vector_general(1.0, A, x, 0.0, y);
  }
} // end PASE_Matrix_multiply_vector()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Matrix_multiply_vector_general"
/**
 * @brief 矩阵向量相乘 y = a * A * x  + b * y (y 不能与 x 相同)
 */
void 
PASE_Matrix_multiply_vector_general(PASE_SCALAR a, PASE_MATRIX A, PASE_VECTOR x, 
                                    PASE_SCALAR b, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Matrix and vectors cannot be empty.\n");
  }
  if((A->global_ncol != x->global_nrow) ||
     (A->global_nrow != y->global_nrow)) {
    PASE_Error(__FUNCT__": The dimensions of matrix and vector are not match.\n");
  }
  if((NULL == A->ops->multiply_matrix_vector_general) &&
     (NULL == A->ops->multiply_matrix_vector)) {
    PASE_Error(__FUNCT__": Neither ops->multiply_matrix_vector_general nor "
                        "ops->multiply_matrix_vector is assigned.\n");
  }
  if(!(A->data_interpreter == x->data_interpreter && A->data_interpreter == y->data_interpreter)) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if(y == x) {
    PASE_Error(__FUNCT__": Vectors x and y cannot be the same.\n");
  }
#endif
  
  if(NULL != A->ops->multiply_matrix_vector_general) {
    A->ops->multiply_matrix_vector_general(a, A->matrix_data, x->vector_data, b, y->vector_data);
  } else {
    PASE_VECTOR tmp_Ax = PASE_Vector_create_by_vector(y);
    PASE_Matrix_multiply_vector(A, x, tmp_Ax);
    PASE_Vector_scale(b, y);
    PASE_Vector_axpy(a, tmp_Ax, y);
    PASE_Vector_destroy(&tmp_Ax);
  }
} // end PASE_Matrix_multiply_vector_general()

#undef  __FUNCT__
#define __FUNCT__ "PASE_MatrixT_multiply_matrix"
/**
 * @brief 矩阵转置与矩阵相乘 A^T * B
 */
PASE_MATRIX 
PASE_MatrixT_multiply_matrix(PASE_MATRIX A, PASE_MATRIX B)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == B)) {
    PASE_Error(__FUNCT__": Neither the tow matrices can be empty.\n");
  }
  if(A->global_nrow != B->global_nrow) {
    PASE_Error(__FUNCT__": Matrix dimensions are not match.\n");
  }
  if(A->data_interpreter != B->data_interpreter) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if((NULL == A->ops->multiply_matrixT_matrix) &&
     ((NULL == A->ops->transpose) || (NULL == A->ops->multiply_matrix_matrix))) {
    PASE_Error(__FUNCT__": Neither ops->multiply_matrixT_matrix nor ops->transpose is assigned.\n");
  }
#endif

  PASE_MATRIX C = NULL;
  if(A->ops->multiply_matrixT_matrix) {
      void *matrix_data = A->ops->multiply_matrixT_matrix(A->matrix_data, B->matrix_data);
      C = PASE_Matrix_assign(matrix_data, A->ops);
      C->is_matrix_data_owner = PASE_YES;
      C->data_interpreter     = A->data_interpreter;

  } else {
      PASE_MATRIX tmp_AT = PASE_Matrix_transpose(A);
      C = PASE_Matrix_multiply_matrix(tmp_AT, B);
      PASE_Matrix_destroy(&tmp_AT);
  }

  return C;
} // end PASE_MatrixT_multiply_matrix()

#undef  __FUNCT__
#define __FUNCT__ "PASE_MatrixT_multiply_vector"
/**
 * @brief 矩阵转置与向量相乘 y = AT * x (其中 y 不能与 x 相同)
 */
void 
PASE_MatrixT_multiply_vector(PASE_MATRIX A, PASE_VECTOR x, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Matrix and vectors cannot be empty.\n");
  }
  if((A->global_nrow != x->global_nrow) ||
     (A->global_ncol != y->global_nrow)) {
    PASE_Error(__FUNCT__": The dimensions of matrix and vector are not matched.\n");
  }
  if((NULL == A->ops->multiply_matrixT_vector) &&
     (NULL == A->ops->multiply_matrixT_vector_general)) {
    PASE_Error(__FUNCT__": Neither ops->multiply_matrixT_vector nor "
                        "ops->multiply_matrixT_vector_general is assigned.\n");
  }
  if(!(A->data_interpreter == x->data_interpreter && A->data_interpreter == y->data_interpreter)) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if(y == x) {
    PASE_Error(__FUNCT__": Vectors x and y cannot be the same.\n");
  }
#endif
  
  if(NULL != A->ops->multiply_matrixT_vector) {
    A->ops->multiply_matrixT_vector(A->matrix_data, x->vector_data, y->vector_data);
  } else {
    PASE_MatrixT_multiply_vector_general(1.0, A, x, 0.0, y);
  }
} // end PASE_MatrixT_multiply_vector()

#undef  __FUNCT__
#define __FUNCT__ "PASE_MatrixT_multiply_vector_general"
/**
 * @brief 矩阵转置与向量相乘 y = a * AT * x  + b * y (其中 y 不能与 x 相同)
 */
void PASE_MatrixT_multiply_vector_general(PASE_SCALAR a, PASE_MATRIX A, PASE_VECTOR x, 
                                          PASE_SCALAR b, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == A) || (NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Matrix and vectors cannot be empty.\n");
  }
  if((A->global_nrow != x->global_nrow) ||
     (A->global_ncol != y->global_nrow)) {
    PASE_Error(__FUNCT__": The dimensions of matrix and vector are not matched.\n");
  }
  if((NULL == A->ops->multiply_matrixT_vector) &&
     (NULL == A->ops->multiply_matrixT_vector_general)) {
    PASE_Error(__FUNCT__": Neither ops->multiply_matrixT_vector nor "
                        "ops->multiply_matrixT_vector_general is assigned.\n");
  }
  if(!(A->data_interpreter == x->data_interpreter && A->data_interpreter == y->data_interpreter)) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if(y == x) {
    PASE_Error(__FUNCT__": Vectors x and y cannot be the same.\n");
  }
#endif
  
  if(NULL != A->ops->multiply_matrixT_vector_general) {
    A->ops->multiply_matrixT_vector_general(1.0, A->matrix_data, x->vector_data, 0.0, y->vector_data);
  } else {
    PASE_VECTOR tmp_ATx = PASE_Vector_create_by_vector(y);
    PASE_MatrixT_multiply_vector(A, x, tmp_ATx);
    PASE_Vector_scale(b, y);
    PASE_Vector_axpy(a, tmp_ATx, y);
    PASE_Vector_destroy(&tmp_ATx);
  }
} // end PASE_MatrixT_multiply_vector_general()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_create_by_matrix_and_vector_data"
/**
 * @brief 由矩阵和向量操作集合创建新的向量
 */
PASE_VECTOR 
PASE_Vector_create_by_matrix_and_vector_data_operator(PASE_MATRIX A, PASE_VECTOR_DATA_OPERATOR ops)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == A) {
    PASE_Error(__FUNCT__": Matrix cannot be empty.\n");
  }
  if(NULL == ops) {
    PASE_Error(__FUNCT__": ops cannot be empty.\n");
  }
  if(NULL == ops->create_by_matrix) {
    PASE_Error(__FUNCT__": ops->create_by_matrix is not assigned.\n");
  }
#endif

  void *vector_data = ops->create_by_matrix(A->matrix_data);
  PASE_VECTOR  x          = PASE_Vector_assign(vector_data, ops);
  x->is_vector_data_owner = PASE_YES;
  x->data_interpreter     = A->data_interpreter;
  
  return x;
} // end PASE_Vector_create_by_matrix_and_vector_data_operator()


#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_set_constant_value"
/**
 * @brief 设置向量元素为常数
 */
void
PASE_Vector_set_constant_value(PASE_VECTOR x, PASE_SCALAR a)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vector cannot be empty.\n");
  }
  if(NULL == x->ops) {
    PASE_Error(__FUNCT__": x->ops is not assigned.\n");
  }
  if(NULL == x->ops->set_constant_value) {
    PASE_Error(__FUNCT__": x->ops->set_constant_value is not assigned.\n");
  }
#endif
  x->ops->set_constant_value(x->vector_data, a);
} // end PASE_Vector_set_constant_value()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_set_random_value"
/**
 * @brief 设置向量元素为随机值
 */
void
PASE_Vector_set_random_value(PASE_VECTOR x, PASE_INT seed)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vector cannot be empty.\n");
  }
  if(NULL == x->ops) {
    PASE_Error(__FUNCT__": x->ops is not assigned.\n");
  }
  if(NULL == x->ops->set_random_value) {
    PASE_Error(__FUNCT__": x->ops->set_random_value is not assigned.\n");
  }
#endif
  x->ops->set_random_value(x->vector_data, seed);
} // end PASE_Vector_set_random_value()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_scale"
/**
 * @brief 向量数乘
 */
void
PASE_Vector_scale(PASE_SCALAR a, PASE_VECTOR x)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vector cannot be empty.\n");
  }
  if(NULL == x->ops) {
    PASE_Error(__FUNCT__": x->ops is not assigned.\n");
  }
  if(NULL == x->ops->scale) {
    PASE_Error(__FUNCT__": x->ops->scale is not assigned.\n");
  }
#endif
  x->ops->scale(a, x->vector_data);
} // end PASE_Vector_scale()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_axpy"
/**
 * @brief 向量校正 y = a * x + y (其中 y 不能与 x 相同)
 */
void
PASE_Vector_axpy(PASE_SCALAR a, PASE_VECTOR x, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Vectors cannot be empty.\n");
  }
  if(NULL == x->ops) {
    PASE_Error(__FUNCT__": x->ops is not assigned.\n");
  }
  if(NULL == x->ops->axpy) {
    PASE_Error(__FUNCT__": x->ops->axpy is not assigned.\n");
  }
  if(x->data_interpreter != y->data_interpreter) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
  if(y == x) {
    PASE_Error(__FUNCT__": Vectors x and y cannot be the same.\n");
  }
#endif
  x->ops->axpy(a, x->vector_data, y->vector_data);
} // end PASE_Vector_axpy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_inner_product"
/**
 * @brief 计算向量内积
 */
PASE_SCALAR
PASE_Vector_get_inner_product(PASE_VECTOR x, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Vectors cannot be empty.\n");
  }
  if(NULL == x->ops) {
    PASE_Error(__FUNCT__": x->ops is not assigned.\n");
  }
  if(NULL == x->ops->get_inner_product) {
    PASE_Error(__FUNCT__": x->ops->get_inner_product is not assigned.\n");
  }
  if(x->data_interpreter != y->data_interpreter) {
    PASE_Error(__FUNCT__": data_interpreters are not match.\n");
  }
#endif

  return x->ops->get_inner_product(x->vector_data, y->vector_data);
} // end PASE_Vector_get_inner_product()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_inner_product_some"
/**
 * @brief 计算向量 x[start], ..., x[end] 全体内积
 */
void
PASE_Vector_get_inner_product_some(PASE_VECTOR *x,
                                   PASE_INT start, PASE_INT end,
                                   PASE_SCALAR **prod)
{
  PASE_INT i = 0;
  PASE_INT j = 0;
  for(i = start; i <= end; ++i) {
    for(j = start; j <= i; ++j) {
      prod[i-start][j-start] = PASE_Vector_get_inner_product(x[j], x[i]);
    }
  }

  for(i = start; i <= end; ++i) {
    for(j = i + 1; j <= end; ++j) {
      prod[i-start][j-start] = prod[j-start][i-start];
    }
  }
} // end PASE_Vector_get_inner_product_some()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_orthogonalize"
/**
 * @brief 将向量 x[i] 与 x[start], ..., x[end] 正交化, i 不能属于 [start, end]
 */
void
PASE_Vector_orthogonalize(PASE_VECTOR *x,
                          PASE_INT i,
                          PASE_INT start, PASE_INT end)
{
#if DEBUG_PASE_ALGEBRA
  if((i >= start) && (i <= end)) {
    PASE_Error(__FUNCT__": index cannot locate in [%d, %d].\n", start, end);
  }
#endif

  PASE_INT  j    = 0;
  PASE_REAL prod = 0.0;
  for(j = start; j <= end; ++j) {
    prod = PASE_Vector_get_inner_product(x[j], x[i]);
    PASE_Vector_axpy(-prod, x[j], x[i]);
  }
  prod = PASE_Vector_get_inner_product(x[i], x[i]);
  assert(!isZero(prod) && "Error @ "__FUNCT__": vector norm cannot be zero.");
  PASE_Vector_scale(1.0/sqrt(prod), x[i]);
} // end PASE_Vector_orthogonalize()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_orthogonalize_all"
/**
 * @brief 将向量 x[0], ..., x[num-1] 全体正交化
 */
void
PASE_Vector_orthogonalize_all(PASE_VECTOR *x, PASE_INT num)
{
  PASE_INT j = 0;
  for(j = 0; j < num; ++j) {
    PASE_Vector_orthogonalize(x, j, 0, j-1);
  }
} // end PASE_Vector_orthogonalize_all()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_combination"
/**
 * @brief 计算向量线性组合: y = \sum_{i=0}^{num_vec-1} coef[i] * x[i]
 *
 * @param x       输入向量, 用于做线性组合的向量组
 * @param num_nev 输入向量, 用于做线型组合的向量个数
 * @param coef    输入向量, 用于做线性组合的系数数组
 * @param y       输出向量, 用于存储线性组合完毕得到的向量
 */
void
PASE_Vector_get_combination(PASE_VECTOR *x, PASE_INT num_vec, 
                            PASE_SCALAR *coef, PASE_VECTOR y)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Cannot compute the linear combination of an empty PASE VECTOR set.\n");
  }
  if(NULL == coef) {
    PASE_Error(__FUNCT__": Cannot compute the linear combination of an empty coefficient set.\n");
  }
  if(NULL == y) {
    PASE_Error(__FUNCT__": Cannot store the linear combination with an empty PASE VECTOR.\n");
  }
#endif

  PASE_INT j = 0;
  PASE_Vector_set_constant_value(y, 0.0);
  for(j = 0; j < num_vec; j++) {
    PASE_Vector_axpy(coef[j], x[j], y);
  }
} // end PASE_Vector_get_combination()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_combination_by_dense_matrix"
/**
 * @brief 计算多重向量线性组合: y = x * mat
 *
 * @param x        输入向量, 用于做线性组合的向量组
 * @param num_nev  输入向量, 用于做线型组合的向量个数
 * @param mat      输入向量, 二维行优先排列数组, 用于做线性组合的系数矩阵, 
 *                          维数为 num_vec * num_mat
 * @param ncol_mat 输入向量, 系数矩阵的列数, 即做线型组合得到的向量个数
 * @param y        输出向量, 用于存储线性组合完毕得到的向量组
 */
void
PASE_Vector_get_combination_by_dense_matrix(PASE_VECTOR *x, PASE_INT num_vec,
                                            PASE_SCALAR **mat, PASE_INT ncol_mat, 
                                            PASE_VECTOR *y)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Cannot compute the linear combinations of an empty PASE VECTOR set.\n");
  }
  if(NULL == mat) {
    PASE_Error(__FUNCT__": Cannot compute the linear combinations of an empty coefficient set.\n");
  }
  if(NULL == y) {
    PASE_Error(__FUNCT__": Cannot store the linear combinations with an empty PASE VECTOR set.\n");
  }
#endif

  PASE_INT i = 0;
  for(i = 0; i < ncol_mat; i++) {
    PASE_Vector_get_combination(x, num_vec, mat[i], y[i]);
  }
} // end PASE_Vector_get_combination_by_dense_matrix()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_inner_product_general"
/**
 * @brief 计算广义向量内积 (x, y)_A
 * @param x  输入向量, 用于计算内积
 * @param y  输入向量, 用于计算内积
 * @param A  输入矩阵, 用于计算内积. 当 A 为 NULL 时退化为普通 l2 内积.
 * @return 广义向量内积
 * @note 不检查矩阵 A 的对称正定性
 */
PASE_SCALAR
PASE_Vector_get_inner_product_general(PASE_VECTOR x, PASE_VECTOR y, PASE_MATRIX A)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == x) || (NULL == y)) {
    PASE_Error(__FUNCT__": Vectors cannot be empty.\n");
  }
#endif

  PASE_SCALAR prod = 0.0;
  if(NULL == A) {
    prod = PASE_Vector_get_inner_product(x, y);
  } else {
    PASE_VECTOR tmp_Ax = PASE_Vector_create_by_vector(x);
    PASE_Matrix_multiply_vector(A, x, tmp_Ax);
    prod = PASE_Vector_get_inner_product(tmp_Ax, y);
    PASE_Vector_destroy(&tmp_Ax);
  }
  return prod;
} // end PASE_Vector_get_inner_product_general()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_get_inner_product_some_general"
/**
 * @brief 计算向量 x[start], ..., x[end] 全体 A 内积
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_get_inner_product_some_general(PASE_VECTOR *x,
                                           PASE_INT start, PASE_INT end, 
                                           PASE_MATRIX A, 
                                           PASE_REAL **prod)
{
#if DEBUG_PASE_ALGEBRA
  if((NULL == x) || (NULL == prod)) {
    PASE_Error(__FUNCT__": Neither vector nor products can be empty.\n");
  }
#endif

  if(NULL == A) {
    PASE_Vector_get_inner_product_some(x, start, end, prod);
  } else {
    PASE_VECTOR tmp_Axi = PASE_Vector_create_by_vector(x[start]);

    PASE_INT i = 0;
    PASE_INT j = 0;
    for(i = start; i <= end; ++i) {
      PASE_Matrix_multiply_vector(A, x[i], tmp_Axi);
      for(j = start; j <= i; ++j) {
        prod[i-start][j-start] = PASE_Vector_get_inner_product(x[j], tmp_Axi);
        prod[j-start][i-start] = prod[i-start][j-start];
      }
    }
    PASE_Vector_destroy(&tmp_Axi);
  }
} // end PASE_Vector_get_inner_product_some_general()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_orthogonalize_general"
/**
 * @brief 将向量 x[i] 与 x[start], ..., x[end] 在 A 内积下正交化, i 不能属于 [start, end]
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_orthogonalize_general(PASE_VECTOR *x, 
                                  PASE_INT i,
                                  PASE_INT start, PASE_INT end, 
                                  PASE_MATRIX A)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vectors cannot be empty.\n");
  }
  if((i >= start) && (i <= end)) {
    PASE_Error(__FUNCT__": index cannot locate in [%d, %d].\n", start, end);
  }
#endif
  
  if(NULL == A) {
    PASE_Vector_orthogonalize(x, i, start, end);
  } else {
    PASE_VECTOR tmp_Axi = PASE_Vector_create_by_vector(x[i]);
    PASE_Matrix_multiply_vector(A, x[i], tmp_Axi);

    PASE_INT  j    = 0;
    PASE_REAL prod = 0.0;
    for(j = start; j <= end; ++j) {
      prod = PASE_Vector_get_inner_product(x[j], tmp_Axi);
      PASE_Vector_axpy(-prod, x[j], x[i]);
    }
    prod = PASE_Vector_get_inner_product(x[i], tmp_Axi);
    assert(!isZero(prod) && "Error @ "__FUNCT__": vector norm cannot be zero.");
    PASE_Vector_scale(1.0/sqrt(prod), x[i]);

    PASE_Vector_destroy(&tmp_Axi);
  }
} // end PASE_Vector_orthogonalize_general()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_orthogonalize_all_general"
/**
 * @brief 将向量 x[0], ..., x[num-1] 全体在 A 内积下正交化
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_orthogonalize_all_general(PASE_VECTOR *x, PASE_INT num, PASE_MATRIX A)
{
#if DEBUG_PASE_ALGEBRA
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vectors cannot be empty.\n");
  }
#endif

  PASE_INT j = 0;
  for(j = 0; j < num; ++j) {
    PASE_Vector_orthogonalize_general(x, j, 0, j-1, A);
  }
} // end PASE_Vector_orthogonalize_all_general()


