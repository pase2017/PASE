#ifndef __PASE_VECTOR_H__
#define __PASE_VECTOR_H__

#include "pase_config.h"

/* 向量操作集合 */
typedef struct PASE_VECTOR_DATA_OPERATOR_PRIVATE_
{
  /**
   * @brief 由向量创建向量
   */
  void * (*create_by_vector)(void *x);

  /**
   * @brief 由矩阵创建向量
   */
  void * (*create_by_matrix)(void *A);

  /**
   * @brief 向量销毁
   */
  void (*destroy)(void *x);

  /**
   * @brief 向量复制
   * @param src  输入参数, 被复制的向量
   * @param dst  输出参数, 复制得到的向量
   */
  void (*copy)(void *src, void *dst);

  /**
   * @brief 设置向量元素为常数
   * @param x  输入/输出参数, 元素被设置为常数的向量
   * @param a  输入参数, 常数
   */
  void (*set_constant_value)(void *x, PASE_SCALAR a);

  /**
   * @brief 设置向量元素为随机数
   * @param x    输入/输出参数, 元素被设置为随机数的向量
   * @param seed 输入参数，随机数种子
   */
  void (*set_random_value)(void *x, PASE_INT seed);

  /**
   * @brief 向量数乘 x = a * x
   * @param a  输入参数,     数乘因子
   * @param x  输入/输出参数, 向量
   */
  void (*scale)(PASE_SCALAR a, void *x);

  /**
   * @brief 向量校正 y = a * x + y
   * @param a  输入参数, 校正参数
   * @param x  输入参数, 校正参数
   * @param y  输入/输出参数, 待校正向量
   */
  void (*axpy)(PASE_SCALAR a, void *x, void *y);

  /**
   * @brief 计算向量内积
   * @param x    输入参数, 待计算内积的向量之一
   * @param y    输入参数, 待计算内积的向量之一
   * @return 向量内积
   * @note 注意到返回值的类型为 PASE_SCALAR, 当涉及复运算时, 需仔细考虑.
   */
  PASE_SCALAR (*get_inner_product)(void *x, void *y);

  /**
   * @brief 获得向量全局长度
   * @param x    输入参数, 待获取全局长度的向量
   * @param nrow 输出参数, 向量全局长度
   */
  PASE_INT (*get_global_nrow)(void *x);

} PASE_VECTOR_DATA_OPERATOR_PRIVATE;

typedef PASE_VECTOR_DATA_OPERATOR_PRIVATE * PASE_VECTOR_DATA_OPERATOR;

typedef struct PASE_VECTOR_PRIVATE_ {
    void                      *vector_data;          // 向量数据
    DATA_INTERPRETER           data_interpreter;     // 矩阵数据解释类型
    PASE_VECTOR_DATA_OPERATOR  ops;                  // 向量运算集合
    PASE_INT                   global_nrow;          // 向量全局长度
    PASE_INT                   is_vector_data_owner; // 向量是否由 PASE 创建
} PASE_VECTOR_PRIVATE;

typedef PASE_VECTOR_PRIVATE * PASE_VECTOR;

/**
 * @brief 指定向量操作集合
 *
 * @param create_by_vector    输入参数, 函数指针, 由向量生成向量
 * @param create_by_matrix    输入参数, 函数指针, 由矩阵生成向量
 * @param destroy             输入参数, 函数指针, 销毁向量数据
 * @param copy                输入参数, 函数指针, 向量数据复制
 * @param set_constant_value  输入参数, 函数指针, 设置向量元素为常数
 * @param set_random_value    输入参数, 函数指针, 设置向量元素为随机数
 * @param scale               输入参数, 函数指针, 向量数乘
 * @param get_inner_product   输入参数, 函数指针, 向量内积
 * @param axpy                输入参数, 函数指针, 向量校正
 * @param get_global_nrow     输入参数, 函数指针, 获取向量全局行数
 *
 * @return PASE_VECTOR_DATA_OPERATOR
 */
PASE_VECTOR_DATA_OPERATOR 
PASE_Vector_data_operator_assign
    (void * (*create_by_vector)  (void *x),
     void * (*create_by_matrix)  (void *A),
     void   (*destroy)           (void *x),
     void   (*copy)              (void *src, void *dst),
     void   (*set_constant_value)(void *x, PASE_SCALAR a),
     void   (*set_random_value)  (void *x, PASE_INT seed),
     void   (*scale)             (PASE_SCALAR a, void *x),
     void   (*axpy)              (PASE_SCALAR a, void *x, void *y),
     PASE_SCALAR (*get_inner_product) (void *x, void *y),
     PASE_INT    (*get_global_nrow)   (void *x));

/**
 * @brief 销毁向量操作集合
 */
void
PASE_Vector_data_operator_destroy(PASE_VECTOR_DATA_OPERATOR *ops);

/**
 * @brief 将向量数据与操作集合关联至向量
 *
 * @param vector_data  输入参数, void * 指针, 指向向量数据
 * @param ops          输入参数, 向量运算集合
 *
 * @return PASE_VECTOR
 */
PASE_VECTOR
PASE_Vector_assign(void *vector_data, PASE_VECTOR_DATA_OPERATOR ops);

/**
 * @brief 销毁向量
 */
void
PASE_Vector_destroy(PASE_VECTOR *x);

/**
 * @brief 复制向量 src 到向量 dst
 */
void
PASE_Vector_copy(PASE_VECTOR src, PASE_VECTOR dst);

/**
 * @brief 获取向量全局行数
 */
PASE_INT
PASE_Vevtor_get_global_nrow(PASE_VECTOR x);

/**
 * @brief 由向量创建向量
 */
PASE_VECTOR
PASE_Vector_create_by_vector(PASE_VECTOR x);

#endif
