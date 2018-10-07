# PASE_MATRIX 结构体

结构体 `PASE_MATRIX` 是 `PASE` 软件的基本数据结构之一，定义矩阵数据及其解释类型、矩阵维数、矩阵操作集合、矩阵创建来源等信息。

``` c
typedef struct PASE_MATRIX_PRIVATE_ {

  void                      *matrix_data;          // 矩阵数据
  DATA_INTERPRETER           data_interpreter;     // 矩阵数据解释类型
  PASE_MATRIX_DATA_OPERATOR  ops;                  // 矩阵操作集合

  PASE_INT                   global_nrow;          // 矩阵全局行数
  PASE_INT                   global_ncol;          // 矩阵全局列数

  PASE_INT                   is_matrix_data_owner; // 矩阵是否由 PASE 创建

} PASE_MATRIX_PRIVATE;

typedef PASE_MATRIX_PRIVATE * PASE_MATRIX;
```

在上述 `PASE_MATRIX` 结构体的定义中，`matrix_data` 为指向实际矩阵数据的 `void *` 指针，由相应矩阵操作集合 `ops` 进行真正的计算；`data_interpreter` 为枚举变量（ `pase_config.h` 定义了 `DATA_INTERPRETER` 枚举类型），用于标记矩阵数据的解释类型（HYPRE, JXPAMG, USER等）。此外，变量 `is_matrix_data_owner` 用于标记矩阵是否由 `PASE` 创建，进而决定释放内存的方式。

矩阵操作集合 `PASE_MATRIX_DATA_OPERATOR` 定义如下。

```c
/* 矩阵运算集合 */
typedef struct PASE_MATRIX_DATA_OPERATOR_PRIVATE_ {
  /**
   * @brief 由矩阵创建矩阵
   */
  void * (*create_by_matrix)(void *A);

  /**
   * @brief 矩阵复制
   * @param src  输入参数, 被复制的矩阵
   * @param dst  输出参数, 复制得到的矩阵
   */
  void (*copy)(void *src, void *dst);

  /**
   * @brief 矩阵销毁
   */
  void (*destroy)(void *A);

  /**
   * @brief 矩阵转置 AT = A
   * @return 转置后的矩阵 AT
   */
  void * (*transpose)(void *A);

  /**
   * @brief 矩阵矩阵相乘 C = A * B
   * @return 矩阵乘积 C
   */
  void * (*multiply_matrix_matrix)  (void *A, void *B);

  /**
   * @brief 矩阵矩阵相乘 C = A^T * B
   * @return 矩阵乘积 C
   */
  void * (*multiply_matrixT_matrix) (void *A, void *B);

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
```

