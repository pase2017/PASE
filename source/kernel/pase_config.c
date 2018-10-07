#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include "pase_config.h"

/*
 * @brief 打印信息
 * @note 仅在 comm 对应的 0 号进程上打印
 */
void PASE_Printf(MPI_Comm comm, char *fmt, ...)
{
  PASE_INT myrank = 0;
  MPI_Comm_rank(comm, &myrank);
  if(0 == myrank) {
    va_list vp;
    va_start(vp, fmt);
    vprintf(fmt, vp);
    va_end(vp);
  }
} // PASE_Printf()

/*
 * @brief 打印错误信息
 */
void PASE_Error(char *fmt, ...)
{
  PASE_INT myrank = 0;
  MPI_Comm_rank(PASE_COMM_WORLD, &myrank);
  if(0 == myrank) {
    printf("PASE ERROR @ ");
    va_list vp;
    va_start(vp, fmt);
    vprintf(fmt, vp);
    va_end(vp);
  }
  MPI_Abort(PASE_COMM_WORLD, -1);
} // PASE_Error()

/*
 * @brief 检测是否为零
 * @note 目前仅支持实型变量
 */
PASE_INT
isZero(PASE_SCALAR a)
{
  if(fabs(a) <= 1.0e-16 ) {
    return 1;
  } else {
    return 0;
  }
} // end isZero()
