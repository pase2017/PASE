#ifndef __PASE_VECTOR_HYPRE_H__
#define __PASE_VECTOR_HYPRE_H__

#include "pase_vector.h"

/**
 * @brief  基于 hypre 生成 PASE 向量操作集合
 * @return PASE_VECTOR_DATA_OPERATOR
 */
PASE_VECTOR_DATA_OPERATOR
PASE_Vector_data_operator_create_by_hypre();

/**
 * @brief  基于 hypre 格式的向量数据生成 PASE 向量
 * @param  vector_data  输入参数, void * 指针, 指向向量数据
 * @return PASE_VECTOR
 */
PASE_VECTOR
PASE_Vector_create_by_hypre(void *vector_data);

#endif
