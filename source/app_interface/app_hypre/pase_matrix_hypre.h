#ifndef __PASE_MATRIX_HYPRE_H__
#define __PASE_MATRIX_HYPRE_H__

#include "pase_matrix.h"

/**
 * @brief  基于 hypre 生成 PASE 矩阵操作集合
 * @return PASE_MATRIX_DATA_OPERATOR
 */
PASE_MATRIX_DATA_OPERATOR
PASE_Matrix_data_operator_create_by_hypre();

/**
 * @brief  基于 hypre 格式的矩阵数据生成 PASE 矩阵
 * @param  matrix_data  输入参数, void * 指针, 指向矩阵数据
 * @return PASE_MATRIX
 */
PASE_MATRIX
PASE_Matrix_create_by_hypre(void *matrix_data);

#endif
