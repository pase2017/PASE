#ifndef __PASE_ALGEBRA_H__
#define __PASE_ALGEBRA_H__
 
#include "pase_config.h"
#include "pase_matrix.h"
#include "pase_vector.h"

/**
 * @brief 矩阵转置
 * @return 转置后的矩阵
 */
PASE_MATRIX
PASE_Matrix_transpose(PASE_MATRIX A);

/**
 * @brief 矩阵与矩阵相乘 A * B
 * @return 矩阵乘积
 */
PASE_MATRIX
PASE_Matrix_multiply_matrix(PASE_MATRIX A, PASE_MATRIX B);

/**
 * @brief 矩阵与向量相乘 y = A * x (其中 y 不能与 x 相同)
 */
void
PASE_Matrix_multiply_vector(PASE_MATRIX A, PASE_VECTOR x, PASE_VECTOR y);

/**
 * @brief 矩阵与向量相乘 y = a * A * x  + b * y (其中 y 不能与 x 相同)
 */
void
PASE_Matrix_multiply_vector_general(PASE_SCALAR a, PASE_MATRIX A, PASE_VECTOR x, 
                                    PASE_SCALAR b, PASE_VECTOR y);

/**
 * @brief 矩阵转置与矩阵相乘 A^T * B
 * @return 矩阵乘积
 */
PASE_MATRIX
PASE_MatrixT_multiply_matrix(PASE_MATRIX A, PASE_MATRIX B);

/**
 * @brief 矩阵转置与向量相乘 y = AT * x (其中 y 不能与 x 相同)
 */
void
PASE_MatrixT_multiply_vector(PASE_MATRIX A, PASE_VECTOR x, PASE_VECTOR y);

/**
 * @brief 矩阵转置与向量相乘 y = a * AT * x  + b * y (其中 y 不能与 x 相同)
 */
void
PASE_MatrixT_multiply_vector_general(PASE_SCALAR a, PASE_MATRIX A, PASE_VECTOR x, 
                                     PASE_SCALAR b, PASE_VECTOR y);

/**
 * @brief 由矩阵和向量操作集合创建新的向量
 * @param A    输入参数, 矩阵
 * @param ops  输入参数, 向量操作集合
 * @return PASE_VECTOR
 */
PASE_VECTOR
PASE_Vector_create_by_matrix_and_vector_data_operator(PASE_MATRIX A, PASE_VECTOR_DATA_OPERATOR ops);

/**
 * @brief 设置向量元素为常数
 */
void
PASE_Vector_set_constant_value(PASE_VECTOR x, PASE_SCALAR a);

/**
 * @brief 设置向量元素为随机值
 */
void
PASE_Vector_set_random_value(PASE_VECTOR x, PASE_INT seed);

/**
 * @brief 向量数乘
 */
void
PASE_Vector_scale(PASE_SCALAR a, PASE_VECTOR x);

/**
 * @brief 向量校正 y = a * x + y (其中 y 不能与 x 相同)
 */
void
PASE_Vector_axpy(PASE_SCALAR a, PASE_VECTOR x, PASE_VECTOR y);

/**
 * @brief 计算向量内积 (x, y)
 */
PASE_SCALAR
PASE_Vector_get_inner_product(PASE_VECTOR x, PASE_VECTOR y);

/**
 * @brief 计算向量 x[start], ..., x[end] 全体内积
 */
void
PASE_Vector_get_inner_product_some(PASE_VECTOR *x,
                                   PASE_INT start, PASE_INT end,
                                   PASE_SCALAR **prod);

/**
 * @brief 将向量 x[i] 与 x[start], ..., x[end] 正交化, i 不能属于 [start, end]
 */
void
PASE_Vector_orthogonalize(PASE_VECTOR *x,
                          PASE_INT i,
                          PASE_INT start, PASE_INT end);

/**
 * @brief 将向量 x[0], ..., x[num-1] 全体正交化
 */
void
PASE_Vector_orthogonalize_all(PASE_VECTOR *x, PASE_INT num);

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
                            PASE_SCALAR *coef, PASE_VECTOR y);

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
                                            PASE_VECTOR *y);

/**
 * @brief 计算广义向量内积 (x, y)_A
 * @param x  输入向量, 用于计算内积
 * @param y  输入向量, 用于计算内积
 * @param A  输入矩阵, 用于计算内积. 当 A 为 NULL 时退化为普通 l2 内积.
 * @return 广义向量内积
 * @note 不检查矩阵 A 的对称正定性
 */
PASE_SCALAR
PASE_Vector_get_inner_product_general(PASE_VECTOR x, PASE_VECTOR y, PASE_MATRIX A);

/**
 * @brief 计算向量 x[start], ..., x[end] 全体 A 内积
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_get_inner_product_some_general(PASE_VECTOR *x,
                                           PASE_INT start, PASE_INT end, 
                                           PASE_MATRIX A, 
                                           PASE_REAL **prod);

/**
 * @brief 将向量 x[i] 与 x[start], ..., x[end] 在 A 内积下正交化, i 不能属于 [start, end]
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_orthogonalize_general(PASE_VECTOR *x, 
                                  PASE_INT i,
                                  PASE_INT start, PASE_INT end, 
                                  PASE_MATRIX A);

/**
 * @brief 将向量 x[0], ..., x[num-1] 全体在 A 内积下正交化
 * @note 不检查矩阵 A 的对称正定性
 */
void
PASE_Vector_orthogonalize_all_general(PASE_VECTOR *x, PASE_INT num, PASE_MATRIX A);

#endif
