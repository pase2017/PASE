#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include "pase_config.h"
#include "pase_matrix_hypre.h"
#include "pase_vector_hypre.h"
#include "pase_algebra.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_ls.h"

#define PRINT_DATA    1
#define NO_PRINT_DATA 0

// get local index
void getLocalIndex(MPI_Comm comm, PASE_INT n, PASE_INT *lower, PASE_INT *upper);
// get the stiffness matrix in 2D Laplacian FDM
void getHypreParCSRMatrix(MPI_Comm comm, HYPRE_IJMatrix *A, PASE_INT n);
// get the unsymmetric modified stiffness matrix in 2D Laplacian FDM
void getHypreParCSRMatrixUnsymmetric(MPI_Comm comm, HYPRE_IJMatrix *A, PASE_INT n);
// print matrix in dense form
void printHypreMatrixData(MPI_Comm comm, HYPRE_ParCSRMatrix A, char *name);
// get vector [1, 1, ..., 1]^T
void getHypreParVectorOne(MPI_Comm comm, HYPRE_IJVector *x, PASE_INT len);
// get vector [1, 2, ..., n]^T
void getHypreParVectorIncrease(MPI_Comm comm, HYPRE_IJVector *x, PASE_INT len);
// print vector in dense form
void printHypreVectorData(MPI_Comm comm, HYPRE_ParVector x, char *name);

void test_getGlobalSize(PASE_MATRIX A);
void test_copyMatrix(PASE_MATRIX A, PASE_INT is_print);
void test_transposeMatrix(PASE_MATRIX A, PASE_INT is_print);
void test_multiplyMatrixMatrix (PASE_MATRIX A, PASE_MATRIX B, PASE_INT is_print);
void test_multiplyMatrixTMatrix(PASE_MATRIX A, PASE_MATRIX B, PASE_INT is_print);

void test_multiplyMatrixVector(PASE_MATRIX A, PASE_VECTOR x, PASE_INT is_print);
void test_multiplyMatrixTVector(PASE_MATRIX A, PASE_VECTOR x, PASE_INT is_print);
void test_orthogonalize_all_general(PASE_MATRIX A, PASE_MATRIX B, PASE_VECTOR x, PASE_INT is_print);

int main(int argc, char *argv[])
{
#if 1
  MPI_Init(&argc, &argv);

  PASE_Printf(MPI_COMM_WORLD, "\n");
  PASE_Printf(MPI_COMM_WORLD, "========================================================\n");
  PASE_Printf(MPI_COMM_WORLD, "======== This is a test program for PASE_MATRIX ========\n");
  PASE_Printf(MPI_COMM_WORLD, "========================================================\n");
  PASE_Printf(MPI_COMM_WORLD, "\n");

  PASE_INT n = 3;
  //GetCommandLineInfo(argc, argv, &n);
  
  // get hypre matrix
  HYPRE_IJMatrix hypreA;
  getHypreParCSRMatrixUnsymmetric(MPI_COMM_WORLD, &hypreA, n);
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJMatrixGetObject(hypreA, (void**) &parcsr_A);

  HYPRE_IJMatrix hypreB;
  getHypreParCSRMatrix(MPI_COMM_WORLD, &hypreB, n);
  HYPRE_ParCSRMatrix parcsr_B;
  HYPRE_IJMatrixGetObject(hypreB, (void**) &parcsr_B);

  PASE_INT num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // get hypre vector
  HYPRE_IJVector hyprex;
  getHypreParVectorOne(MPI_COMM_WORLD, &hyprex, n*n);
  HYPRE_ParVector par_x;
  HYPRE_IJVectorGetObject(hyprex, (void **) &par_x);

  HYPRE_IJVector hyprey;
  getHypreParVectorIncrease(MPI_COMM_WORLD, &hyprey, n*n);
  HYPRE_ParVector par_y;
  HYPRE_IJVectorGetObject(hyprey, (void **) &par_y);

  // get PASE matrix
  PASE_MATRIX A = PASE_Matrix_create_by_hypre((void*)parcsr_A);
  PASE_MATRIX B = PASE_Matrix_create_by_hypre((void*)parcsr_B);

  // get PASE vector
  PASE_VECTOR x = PASE_Vector_create_by_hypre((void*)par_x);
  PASE_VECTOR y = PASE_Vector_create_by_hypre((void*)par_y);

  /*---------- test begin ----------*/
  test_getGlobalSize(A);
  test_copyMatrix(A, PRINT_DATA);
  test_transposeMatrix(A, PRINT_DATA);
  test_multiplyMatrixMatrix (A, B, PRINT_DATA);
  test_multiplyMatrixTMatrix(A, B, PRINT_DATA);

  test_multiplyMatrixVector(A, x, PRINT_DATA);
  test_multiplyMatrixTVector(A, x, PRINT_DATA);
  test_orthogonalize_all_general(A, B, x, PRINT_DATA);
  test_orthogonalize_all_general(A, NULL, x, PRINT_DATA);
  /*---------- test end ----------*/

  PASE_Vector_destroy(&y);
  PASE_Vector_destroy(&x);
  PASE_Matrix_destroy(&B);
  PASE_Matrix_destroy(&A);

  HYPRE_IJVectorDestroy(hyprey);
  HYPRE_IJVectorDestroy(hyprex);
  HYPRE_IJMatrixDestroy(hypreB);
  HYPRE_IJMatrixDestroy(hypreA);

  MPI_Finalize();
#endif
  return 0;
}

void
getLocalIndex(MPI_Comm comm, PASE_INT n, PASE_INT *lower, PASE_INT *upper)
{
  PASE_INT myid, num_procs;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &num_procs);

  /* Each processor knows only of its own rows - the range is denoted by 
   * lower and upper.  We account for the fact that the row number 
   * may not divide evenly by the number of processors. 
   **/
  PASE_INT local_size, extra;

  local_size = n / num_procs;
  extra      = n - local_size * num_procs;

  *lower  = local_size * myid;
  *lower += hypre_min(myid, extra);

  *upper  = local_size * (myid+1);
  *upper += hypre_min(myid+1, extra);
  *upper  = *upper - 1;
}

void
getHypreParCSRMatrix(MPI_Comm comm, HYPRE_IJMatrix *A, PASE_INT n)
{
  PASE_INT dim = n * n;
  PASE_INT ilower = -1;
  PASE_INT iupper = -1;
  getLocalIndex(comm, dim, &ilower, &iupper);

  /* Create matrix.
   * Note that this is a square matrix, so we indicate the row partition
   * size twice (since number of rows = number of cols).
   **/
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, A);
  /* Choose a parallel csr format storage */
  HYPRE_IJMatrixSetObjectType(*A, HYPRE_PARCSR);
  /* Initialize before setting coefficients */
  HYPRE_IJMatrixInitialize(*A);

  /* Now go through my local rows and set the matrix entries.
   * Each row has at most 5 entries. For example, if n=3:
   *   A = [M -I 0; -I M -I; 0 -I M]
   *   M = [4 -1 0; -1 4 -1; 0 -1 4]
   **/
  PASE_INT  nnz;
  PASE_REAL values[5];
  PASE_INT  cols[5];
  PASE_INT  i;
  for(i = ilower; i <= iupper; i++) {
    nnz = 0;

    /* The left -1 in -I : position i-n */
    if( i >= n)
    {
      cols[nnz]   = i - n;
      values[nnz] = -1.0;
      nnz++;
    }

    /* The left -1 in M : position i-1 */
    if(i % n)
    {
      cols[nnz]   = i - 1;
      values[nnz] = -1.0;
      nnz++;
    }

    /* The diagonal : position i */
    cols[nnz]   = i;
    values[nnz] = 4.0;
    nnz++;

    /* The right -1 in M : position i+1 */
    if((i+1) % n)
    {
      cols[nnz]  = i + 1;
      values[nnz] = -1.0;
      nnz++;
    }

    /* The right -1 in -I : position i+n */
    if((i+n) < dim)
    {
      cols[nnz]   = i + n;
      values[nnz] = -1.0;
      nnz++;
    }

    /* Set the values for row i */
    HYPRE_IJMatrixSetValues(*A, 1, &nnz, &i, cols, values);
  }

  HYPRE_IJMatrixAssemble(*A);
} // end getHypreParCSRMatrixStiff()

void
getHypreParCSRMatrixUnsymmetric(MPI_Comm comm, HYPRE_IJMatrix *A, PASE_INT n)
{
  PASE_INT dim = n * n;
  PASE_INT ilower = -1;
  PASE_INT iupper = -1;
  getLocalIndex(comm, dim, &ilower, &iupper);

  /* Create matrix.
   * Note that this is a square matrix, so we indicate the row partition
   * size twice (since number of rows = number of cols).
   **/
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, A);
  /* Choose a parallel csr format storage */
  HYPRE_IJMatrixSetObjectType(*A, HYPRE_PARCSR);
  /* Initialize before setting coefficients */
  HYPRE_IJMatrixInitialize(*A);

  /* Now go through my local rows and set the matrix entries.
   * Each row has at most 5 entries. For example, if n=3:
   *   A = [M -3I 0; -I M -3I; 0 -I M]
   *   M = [4 -3  0; -1 4 -3;  0 -1 4]
   **/
  PASE_INT  nnz;
  PASE_REAL values[5];
  PASE_INT  cols[5];
  PASE_INT  i;
  for(i = ilower; i <= iupper; i++) {
    nnz = 0;

    /* The left -1 in -I : position i-n */
    if( i >= n)
    {
      cols[nnz]   = i - n;
      values[nnz] = -1.0;
      nnz++;
    }

    /* The left -1 in M : position i-1 */
    if(i % n)
    {
      cols[nnz]   = i - 1;
      values[nnz] = -1.0;
      nnz++;
    }

    /* The diagonal : position i */
    cols[nnz]   = i;
    values[nnz] = 4.0;
    nnz++;

    /* The right -3 in M : position i+1 */
    if((i+1) % n)
    {
      cols[nnz]  = i + 1;
      values[nnz] = -3.0;
      nnz++;
    }

    /* The right -3 in -3I : position i+n */
    if((i+n) < dim)
    {
      cols[nnz]   = i + n;
      values[nnz] = -3.0;
      nnz++;
    }

    /* Set the values for row i */
    HYPRE_IJMatrixSetValues(*A, 1, &nnz, &i, cols, values);
  }

  HYPRE_IJMatrixAssemble(*A);
} // end getHypreParCSRMatrixStiffLower()

void printHypreMatrixData(MPI_Comm comm, HYPRE_ParCSRMatrix A, char *name)
{
  PASE_Printf(comm, "\n  %s\n", name);

  PASE_INT myid, num_procs;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &num_procs);

  int num_local_row = A->last_row_index - A->first_row_index + 1;
  int idx_proc = 0;
  for(idx_proc = 0; idx_proc < num_procs; ++idx_proc) {
    if(myid == idx_proc) {
      int idx_local_row = 0;
      for(idx_local_row = 0; idx_local_row < num_local_row; ++idx_local_row) {
        printf("    proc % 2d:  ", myid);
        int idx_col = 0;
        for(idx_col = 0; idx_col < A->global_num_cols; ++idx_col) {
          int is_col_non_zero = 0;
          int row_start_diag  = A->diag->i[idx_local_row];
          int row_end_diag    = A->diag->i[idx_local_row+1];
          int row_start_offd  = A->offd->i[idx_local_row];
          int row_end_offd    = A->offd->i[idx_local_row+1];
          int idx_local_col   = 0;
          for(idx_local_col = row_start_diag; idx_local_col < row_end_diag; ++idx_local_col) {
            if(idx_col == (A->diag->j[idx_local_col] + A->first_row_index)) {
              printf("% 6.2f ", A->diag->data[idx_local_col]);
              is_col_non_zero = 1;
              break;
            }
          }
          for(idx_local_col = row_start_offd; idx_local_col<row_end_offd; ++idx_local_col) {
            if(idx_col == A->col_map_offd[A->offd->j[idx_local_col]]) {
              printf("% 6.2f ", A->offd->data[idx_local_col]);
              is_col_non_zero = 1;
              break;
            }
          }
          if(0 == is_col_non_zero) {
            printf("% 6.2f ", 0.0);
          }
        }
        printf("\n");
        fflush(stdout);
      } 
    } else {
      usleep(100);
    }
    MPI_Barrier(comm);
  }
  PASE_Printf(comm, "\n");
} // end printHypreMatrixData()

// get hypre vector [1, 1, ..., 1]^T
void getHypreParVectorOne(MPI_Comm comm, HYPRE_IJVector *x, PASE_INT len)
{
  PASE_INT ilower = -1;
  PASE_INT iupper = -1;
  getLocalIndex(comm, len, &ilower, &iupper);
  PASE_INT local_size   = iupper - ilower + 1;

  HYPRE_IJVectorCreate(comm, ilower, iupper, x);
  HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*x);
 
  PASE_SCALAR *x_values = (PASE_SCALAR*)calloc(local_size, sizeof(PASE_SCALAR));
  PASE_INT    *rows     = (PASE_INT*)   calloc(local_size, sizeof(PASE_INT));

  PASE_INT i = 0;
  for(i=0; i<local_size; i++) {
    x_values[i] = 1.0;
    rows[i]     = ilower + i;
  }
  
  HYPRE_IJVectorSetValues(*x, local_size, rows, x_values);
  
  free(x_values);
  free(rows);

  HYPRE_IJVectorAssemble(*x);
} // end getHypreParVectorOne()

// get hypre vector [1, 2, ..., n]^T
void getHypreParVectorIncrease(MPI_Comm comm, HYPRE_IJVector *x, PASE_INT len)
{
  PASE_INT ilower = -1;
  PASE_INT iupper = -1;
  getLocalIndex(comm, len, &ilower, &iupper);
  PASE_INT local_size   = iupper - ilower + 1;

  HYPRE_IJVectorCreate(comm, ilower, iupper,x);
  HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*x);
 
  PASE_SCALAR *x_values = (PASE_SCALAR*)calloc(local_size, sizeof(PASE_SCALAR));
  PASE_INT    *rows     = (PASE_INT*)   calloc(local_size, sizeof(PASE_INT));

  PASE_INT i = 0;
  for(i=0; i<local_size; i++) {
    x_values[i] = ilower + i;
    rows[i]     = ilower + i;
  }
  
  HYPRE_IJVectorSetValues(*x, local_size, rows, x_values);
  
  free(x_values);
  free(rows);

  HYPRE_IJVectorAssemble(*x);
} // end getHypreParVectorIncrease()

void printHypreVectorData(MPI_Comm comm, HYPRE_ParVector x, char *name)
{
  PASE_Printf(comm, "\n  %s\n", name);

  PASE_INT myid, num_procs;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &num_procs);

  PASE_INT local_size = x->local_vector->size;
  PASE_INT idx_proc = 0;
  for(idx_proc = 0; idx_proc < num_procs; idx_proc++) {
    if(myid == idx_proc) {
      PASE_INT i = 0;
      for(i = 0; i<local_size; i++) {
        printf("    proc % 2d: % 6.2f \n", myid, x->local_vector->data[i]);
      }
      fflush(stdout);
    } else {
	    usleep(100);
    }
    MPI_Barrier(comm);
  }
  PASE_Printf(comm, "\n");
} // end printHypreVectorData()

#undef  __FUNCT__
#define __FUNCT__ "test_getGlobalSize"
void test_getGlobalSize(PASE_MATRIX A)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  PASE_INT nrow = PASE_Matrix_get_global_nrow(A);
  PASE_INT ncol = PASE_Matrix_get_global_ncol(A);
  PASE_Printf(comm, "    Matrix size: %d x %d\n\n", nrow, ncol);
} // end test_getGlobalSize()

#undef  __FUNCT__
#define __FUNCT__ "test_copyMatrix"
void test_copyMatrix(PASE_MATRIX A, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");

  int nrow = PASE_Matrix_get_global_nrow(A);
  HYPRE_IJMatrix hypreB;
  getHypreParCSRMatrix(comm, &hypreB, sqrt(nrow));
  HYPRE_ParCSRMatrix parcsr_B;
  HYPRE_IJMatrixGetObject(hypreB, (void**) &parcsr_B);
  PASE_MATRIX B = PASE_Matrix_create_by_hypre((void*)parcsr_B);
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)B->matrix_data, "PASE B");

  PASE_Matrix_copy(A, B);
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)B->matrix_data, "PASE B copied");

  PASE_Matrix_destroy(&B);
  HYPRE_IJMatrixDestroy(hypreB);
} // end test_copyMatrix()

#undef  __FUNCT__
#define __FUNCT__ "test_transposeMatrix"
void test_transposeMatrix(PASE_MATRIX A, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");
  PASE_MATRIX AT = PASE_Matrix_transpose(A);
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)AT->matrix_data, "PASE AT");
  PASE_Matrix_destroy(&AT);
} // end test_transposeMatrix()

#undef  __FUNCT__
#define __FUNCT__ "test_multiplyMatrixMatrix"
void test_multiplyMatrixMatrix (PASE_MATRIX A, PASE_MATRIX B, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)B->matrix_data, "PASE B");
  PASE_MATRIX AB = PASE_Matrix_multiply_matrix(A, B);
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)AB->matrix_data, "PASE AB");
  PASE_Matrix_destroy(&AB);
} // end test_multiplyMatrixMatrix()

#undef  __FUNCT__
#define __FUNCT__ "test_multiplyMatrixTMatrix"
void test_multiplyMatrixTMatrix(PASE_MATRIX A, PASE_MATRIX B, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)B->matrix_data, "PASE B");
  PASE_MATRIX ATB = PASE_MatrixT_multiply_matrix(A, B);
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)ATB->matrix_data, "PASE ATB");
  PASE_Matrix_destroy(&ATB);
} // end test_multiplyMatrixTMatrix()

#undef  __FUNCT__
#define __FUNCT__ "test_multiplyMatrixVector"
void test_multiplyMatrixVector(PASE_MATRIX A, PASE_VECTOR x, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)   x->vector_data, "PASE x");

  PASE_VECTOR y = PASE_Vector_create_by_matrix_and_vector_data_operator(A, x->ops);

  PASE_Matrix_multiply_vector(A, x, y);
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)y->vector_data, "PASE y = A * x");

  PASE_Matrix_multiply_vector_general(1.0, A, x, 0.0, y);
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)y->vector_data, "PASE y = 1.0 * A * x + 0.0 * y");

  PASE_Vector_destroy(&y);
} // end test_multiplyMatrixVector()

#undef  __FUNCT__
#define __FUNCT__ "test_multiplyMatrixTVector"
void test_multiplyMatrixTVector(PASE_MATRIX A, PASE_VECTOR x, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");
  if(is_print) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)A->matrix_data, "PASE A");
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)   x->vector_data, "PASE x");

  PASE_VECTOR y = PASE_Vector_create_by_matrix_and_vector_data_operator(A, x->ops);

  PASE_MatrixT_multiply_vector(A, x, y);
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)y->vector_data, "PASE y = AT * x");

  PASE_MatrixT_multiply_vector_general(1.0, A, x, 0.0, y);
  if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)y->vector_data, "PASE y = 1.0 * AT * x + 0.0 * y");

  PASE_Vector_destroy(&y);
} // end test_multiplyMatrixTVector()

#undef  __FUNCT__
#define __FUNCT__ "test_orthogonalize_all_general"
void test_orthogonalize_all_general(PASE_MATRIX A, PASE_MATRIX B, PASE_VECTOR x, PASE_INT is_print)
{
  MPI_Comm comm = PASE_Matrix_get_mpi_comm(A);
  PASE_Printf(comm, "## "__FUNCT__"\n");

  if(is_print && NULL != B) printHypreMatrixData(comm, (HYPRE_ParCSRMatrix)B->matrix_data, "PASE B");

  PASE_INT len = 4;
  PASE_VECTOR *vec = (PASE_VECTOR*)PASE_Malloc(len * sizeof(PASE_VECTOR));
  PASE_INT i = 0, j = 0;
  vec[0] = PASE_Vector_create_by_vector(x);
  PASE_Vector_set_constant_value(vec[0], 2.0);
  //if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)vec[0]->vector_data, "PASE vec[0]");
  for(i = 1; i < len; ++i) {
    vec[i] = PASE_Vector_create_by_matrix_and_vector_data_operator(A, x->ops);
    PASE_Vector_set_random_value(vec[i], i);
    //if(is_print) printHypreVectorData(comm, (HYPRE_ParVector)vec[i]->vector_data, "PASE vec[i]");
  }

  PASE_REAL **prod = (PASE_REAL**)PASE_Malloc(len * sizeof(PASE_REAL*));
  for(i = 0; i < len; ++i) {
    prod[i] = (PASE_REAL*)PASE_Malloc(len * sizeof(PASE_REAL));
  }
  if(is_print) {
    PASE_Vector_get_inner_product_some_general(vec, 0, len-1, B, prod);
    for(i = 0; i < len; ++i) {
      PASE_Printf(comm, "    prod before ortho:  ");
      for(j = 0; j < len; ++j) {
        PASE_Printf(comm, " % 6.2f ", prod[i][j]);
      }
      PASE_Printf(comm, "\n");
    }
    PASE_Printf(comm, "\n");
  }

  PASE_Vector_orthogonalize_all_general(vec, len, B);
  if(is_print) {
    PASE_Vector_get_inner_product_some_general(vec, 0, len-1, B, prod);
    for(i = 0; i < len; ++i) {
      PASE_Printf(comm, "    prod after  ortho:  ");
      for(j = 0; j < len; ++j) {
        PASE_Printf(comm, " % 6.2f ", prod[i][j]);
      }
      PASE_Printf(comm, "\n");
    }
    PASE_Printf(comm, "\n");
  }
  
  for(i = 0; i < len; ++i) {
    PASE_Free(prod[i]);
  }
  PASE_Free(prod);
  for(i = 0; i < len; ++i) {
    PASE_Vector_destroy(vec + i);
  }
  PASE_Free(vec);
} // end test_orthogonalize_all_general()