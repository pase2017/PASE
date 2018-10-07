#ifndef __PASE_MATRIX_H__
#define __PASE_MATRIX_H__
 
#include "pase_config.h"

/* 矩阵操作集合 */
typedef struct PASE_MATRIX_DATA_OPERATOR_PRIVATE_ 
{
  /**
   * @brief 由矩阵创建矩阵
   */
  void * (*create_by_matrix)(void *A);

  /**
   * @brief 矩阵销毁
   */
  void (*destroy)(void *A);

  /**
   * @brief 矩阵复制
   * @param src  输入参数, 被复制的矩阵
   * @param dst  输出参数, 复制得到的矩阵
   * @note 1. 该函数不负责分配内存, 因此 src 和 dst 的 matrix_data 不能为空
   *       2. 假设 src 和 dst 具有相同的稀疏结构
   */
  void (*copy)(void *src, void *dst);

  /**
   * @brief 矩阵转置
   * @return 转置后的矩阵
   */
  void * (*transpose)(void *A);

  /**
   * @brief 矩阵矩阵相乘 C = A * B
   * @return 矩阵乘积 C
   */
  void * (*multiply_matrix_matrix)  (void *A, void *B);

  /**
   * @brief 矩阵向量相乘 y = A * x 
   */
  void (*multiply_matrix_vector)(void *A, void *x, void *y);

  /**
   * @brief 矩阵向量相乘 y = a * A * x  + b * y
   */
  void (*multiply_matrix_vector_general)(PASE_SCALAR a, void *A, void *x, 
                                         PASE_SCALAR b, void *y);

  /**
   * @brief 矩阵矩阵相乘 C = A^T * B
   * @return 矩阵乘积 C
   */
  void * (*multiply_matrixT_matrix) (void *A, void *B);

  /**
   * @brief 矩阵向量相乘 y = A^T * x 
   */
  void (*multiply_matrixT_vector) (void *A, void *x, void *y);

  /**
   * @brief 矩阵向量相乘 y = a * A^T * x  + b * y
   */
  void (*multiply_matrixT_vector_general)(PASE_SCALAR a, void *A, void *x, 
                                          PASE_SCALAR b, void *y);

  /**
   * @brief 获得矩阵全局行数
   */
  PASE_INT (*get_global_nrow)(void *A);

  /**
   * @brief 获得矩阵全局列数
   */
  PASE_INT (*get_global_ncol)(void *A);

  /**
   * @brief 获得矩阵通信器
   */
  MPI_Comm (*get_mpi_comm)(void *A);

} PASE_MATRIX_DATA_OPERATOR_PRIVATE;

typedef PASE_MATRIX_DATA_OPERATOR_PRIVATE * PASE_MATRIX_DATA_OPERATOR;

// PASE 矩阵结构
typedef struct PASE_MATRIX_PRIVATE_ 
{
  void                      *matrix_data;          // 矩阵数据
  DATA_INTERPRETER           data_interpreter;     // 矩阵数据解释类型
  PASE_MATRIX_DATA_OPERATOR  ops;                  // 矩阵操作集合

  PASE_INT                   global_nrow;          // 矩阵全局行数
  PASE_INT                   global_ncol;          // 矩阵全局列数

  PASE_INT                   is_matrix_data_owner; // 矩阵是否由 PASE 创建

} PASE_MATRIX_PRIVATE;

typedef PASE_MATRIX_PRIVATE * PASE_MATRIX;

/**
 * @brief 指定矩阵操作集合
 *      
 * @param create_by_matrix                 输入参数, 函数指针, 由矩阵生成矩阵
 * @param destroy                          输入参数, 函数指针, 销毁矩阵数据
 * @param copy                             输入参数, 函数指针, 矩阵复制
 * @param transpose                        输入参数, 函数指针, 矩阵转置
 * @param multiply_matrix_matrix           输入参数, 函数指针, 矩阵与矩阵相乘
 * @param multiply_matrix_vector           输入参数, 函数指针, 矩阵与向量相乘
 * @param multiply_matrix_vector_general   输入参数, 函数指针, 矩阵与向量相乘并与另一向量做线性组合
 * @param multiply_matrixT_matrix          输入参数, 函数指针, 矩阵转置与矩阵相乘
 * @param multiply_matrixT_vector          输入参数, 函数指针, 矩阵转置与向量相乘
 * @param multiply_matrixT_vector_general  输入参数, 函数指针, 矩阵转置与向量相乘并与另一向量做线性组合
 * @param get_global_nrow                  输入参数, 函数指针, 获取矩阵全局行数
 * @param get_global_ncol                  输入参数, 函数指针, 获取矩阵全局列数
 * @param get_mpi_comm                     输入参数, 函数指针, 获取矩阵通信器
 *
 * @note  multiply_matrix_vector  与 multiply_matrix_vector_general 不能同时为 NULL
 * @note  multiply_matrixT_vector 与 multiply_matrixT_vector_general 不能同时为 NULL
 *
 * @return PASE_VECTOR_DATA_OPERATOR
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
     MPI_Comm (*get_mpi_comm)                 (void *A));

/**
 * @brief 销毁矩阵操作集合
 */
void
PASE_Matrix_data_operator_destroy(PASE_MATRIX_DATA_OPERATOR *ops);

/**
 * @brief 将矩阵数据与操作集合关联至矩阵
 *
 * @param matrix_data  输入参数, void * 指针, 指向矩阵数据
 * @param ops          输入参数, 矩阵操作集合
 *
 * @return PASE_MATRIX
 */
PASE_MATRIX
PASE_Matrix_assign(void *matrix_data, PASE_MATRIX_DATA_OPERATOR ops);

/**
 * @brief 销毁矩阵
 */
void
PASE_Matrix_destroy(PASE_MATRIX *A);

/**
 * @brief 矩阵复制
 * @param src  输入参数, 被复制的矩阵
 * @param dst  输出参数, 复制得到的矩阵
 * @note 1. 该函数不负责分配内存, 因此 src 和 dst 的 matrix_data 不能为空
 *       2. 仅用于 src 和 dst 具有相同稀疏结构的情况
 */
void
PASE_Matrix_copy(PASE_MATRIX src, PASE_MATRIX dst);

/**
 * @brief 获取矩阵全局行数
 */
PASE_INT
PASE_Matrix_get_global_nrow(PASE_MATRIX A);

/**
 * @brief 获取矩阵全局列数
 */
PASE_INT
PASE_Matrix_get_global_ncol(PASE_MATRIX A);

/**
 * @brief 获取矩阵通信器
 */
MPI_Comm
PASE_Matrix_get_mpi_comm(PASE_MATRIX A);

#endif
