#ifndef __PASE_CONFIG_H__
#define __PASE_CONFIG_H__

//=============================================================================
/* 基本数据类型 */
typedef int    PASE_INT;
typedef double PASE_DOUBLE;
typedef double PASE_REAL;

#ifdef PASE_USE_COMPLEX
typedef complex PASE_SCALAR;
typedef complex PASE_COMPLEX;
#else
typedef double  PASE_SCALAR;
#endif

//=============================================================================
/* 是否使用 HYPRE 软件包 */
#ifdef USE_HYPRE
#define PASE_USE_HYPRE 1 
#else
#define PASE_USE_HYPRE 0 
#endif

/* 是否使用 JXPAMG 软件包 */
#ifdef USE_JXPAMG
#define PASE_USE_JXPAMG 1
#else
#define PASE_USE_JXPAMG 0
#endif

//=============================================================================
/* 预定义逻辑变量 */
#define PASE_NO       0
#define PASE_YES      1
#define PASE_UNKNOWN -1

//=============================================================================
/* 预定义枚举类型 */
typedef enum { CLJP = 1, FALGOUT = 2, PMHIS = 3 } PASE_COARSEN_TYPE;
typedef enum { DATA_INTERPRETER_USER   = -1, 
               DATA_INTERPRETER_HYPRE  =  1, 
               DATA_INTERPRETER_JXPAMG =  2 } DATA_INTERPRETER;

//=============================================================================
/* MPI 以及打印函数 */
#include "mpi.h"

#define PASE_COMM_WORLD MPI_COMM_WORLD
#define PASE_COMM_SELF  MPI_COMM_SELF

/*
 * @brief 打印错误信息
 */
void PASE_Error(char *fmt, ...);

/*
 * @brief 打印信息
 * @note 仅在 comm 对应的 0 号进程上打印
 */
void PASE_Printf(MPI_Comm comm, char *fmt, ...);

//=============================================================================
/* 内存管理 */
#include <stdlib.h>

#define PASE_Malloc  malloc
#define PASE_Realloc realloc
#define PASE_Calloc  calloc
#define PASE_Free(a) { free(a); a = NULL; }

//=============================================================================
/*
 * @brief 检测是否为零
 * @note 目前仅支持实型变量
 */
PASE_INT
isZero(PASE_SCALAR a);

#endif
