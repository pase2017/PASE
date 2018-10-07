#include "pase_config.h"
#include "pase_vector.h"

#define DEBUG_PASE_VECTOR 1

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_data_operator_assign"
/**
 * @brief 指定向量操作集合
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
     PASE_INT    (*get_global_nrow)   (void *x))
{
  PASE_VECTOR_DATA_OPERATOR ops
      = (PASE_VECTOR_DATA_OPERATOR)PASE_Malloc(sizeof(PASE_VECTOR_DATA_OPERATOR_PRIVATE));

  ops->create_by_vector   = create_by_vector;
  ops->create_by_matrix   = create_by_matrix;
  ops->destroy            = destroy;
  ops->copy               = copy;
  ops->set_constant_value = set_constant_value;
  ops->set_random_value   = set_random_value;
  ops->scale              = scale;
  ops->axpy               = axpy;
  ops->get_inner_product  = get_inner_product;
  ops->get_global_nrow    = get_global_nrow;

  return ops;
} // end PASE_Vector_data_operator_assign()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_data_operator_destroy"
/**
 * @brief 销毁向量操作集合
 */
void
PASE_Vector_data_operator_destroy(PASE_VECTOR_DATA_OPERATOR *ops)
{
  PASE_Free(*ops);
} // end PASE_Vector_data_operator_destroy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_assign"
/**
 * @brief 将向量数据与操作集合关联至向量
 */
PASE_VECTOR
PASE_Vector_assign(void *vector_data, PASE_VECTOR_DATA_OPERATOR ops)
{
#if DEBUG_PASE_VECTOR
  if(NULL == vector_data) {
    PASE_Error(__FUNCT__": Cannot create PAES VECTOR without vector data.\n");
  }
  if(NULL == ops) {
    PASE_Error(__FUNCT__": Cannot create PAES VECTOR without vector data operator.\n");
  }
#endif

  PASE_VECTOR x = (PASE_VECTOR)PASE_Malloc(sizeof(PASE_VECTOR_PRIVATE));
  x->ops = (PASE_VECTOR_DATA_OPERATOR)PASE_Malloc(sizeof(PASE_VECTOR_DATA_OPERATOR_PRIVATE));
  *(x->ops)               = *ops;
  x->vector_data          = vector_data;
  x->is_vector_data_owner = PASE_NO;
  x->data_interpreter     = DATA_INTERPRETER_USER;
  x->global_nrow          = x->ops->get_global_nrow(x->vector_data);
  return x;
} // end PASE_Vector_assign()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_destroy"
/**
 * @brief 销毁向量
 */
void
PASE_Vector_destroy(PASE_VECTOR *x)
{
  if(NULL == *x) return;

#if DEBUG_PASE_VECTOR
  if((PASE_YES != (*x)->is_vector_data_owner) &&
     (PASE_NO  != (*x)->is_vector_data_owner)) {
    PASE_Error(__FUNCT__": Cannot decide whether the owner of vector is.");
  }
  if(PASE_YES == (*x)->is_vector_data_owner) {
    if(NULL == (*x)->ops->destroy) {
      PASE_Error(__FUNCT__": Operator ops->destroy is not assigned.\n");
    }
  }
#endif

  if(PASE_YES == (*x)->is_vector_data_owner) {
    (*x)->ops->destroy((*x)->vector_data);
  }
  PASE_Vector_data_operator_destroy(&((*x)->ops));
  PASE_Free(*x);
} // end PASE_Vector_destroy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_copy"
/**
 * @brief 复制向量 src 到向量 dst
 */
void
PASE_Vector_copy(PASE_VECTOR src, PASE_VECTOR dst)
{
#if DEBUG_PASE_VECTOR
  if(src->global_nrow != dst->global_nrow) {
    PASE_Error(__FUNCT__": Cannot copy vector if dimensions are not matched.\n");
  }
#endif

  src->ops->copy(src->vector_data, dst->vector_data);
} // end PASE_Vector_copy()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vevtor_get_global_nrow"
/**
 * @brief 获取向量全局行数
 */
PASE_INT
PASE_Vevtor_get_global_nrow(PASE_VECTOR x)
{
#if DEBUG_PASE_VECTOR
  if(NULL == x) {
    PASE_Error(__FUNCT__": Vector cannot be empty.\n");
  }
  if(NULL == x->ops->get_global_nrow) {
    PASE_Error(__FUNCT__": Operator ops->get_global_nrow is not assigned.\n");
  }
#endif
  return x->ops->get_global_nrow(x->vector_data);
} // end PASE_Vevtor_get_global_nrow()

#undef  __FUNCT__
#define __FUNCT__ "PASE_Vector_create_by_vector"
/**
 * @brief 由向量创建向量
 */
PASE_VECTOR
PASE_Vector_create_by_vector(PASE_VECTOR x)
{
  void        *vector_data = x->ops->create_by_vector(x->vector_data);
  PASE_VECTOR  y           = PASE_Vector_assign(vector_data, x->ops);
  y->is_vector_data_owner  = PASE_YES;
  y->data_interpreter      = x->data_interpreter;
  return y;
} // end PASE_Vector_create_by_vector()