/**
 * @file bblas_zmkl.c
 *
 * @brief BBLAS MKL routines for double _Complex precision.
 *
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-03-17
 *
 * Contains wrappers to batch MKL functions.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include "bblas_mkl.h"

/**
 * Wrapper for zgemm_batch in MKL.
 *
 * Calls cblas_zgemm_batch (linked to MKL) with the input parameters and prints out
 * any errors that may occur.
 *
 * The calling sequence is identical to our reference implementation of BBLAS.
 * The arguments are converted to an MKL function call within this function.
 **/
void mkl_zgemm_batch(const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
					 const int *M,  const int *N, const int *K,
					 const BBLAS_Complex64_t *alpha,
					 const BBLAS_Complex64_t **arrayA, const int *lda,
					 const BBLAS_Complex64_t **arrayB,
					 const int *ldb, const BBLAS_Complex64_t *beta,
					 BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
					 enum BBLAS_OPTS batch_opts, int *info)
{
#if defined(BBLAS_WITH_MKL)

	if (batch_opts == BBLAS_VARIABLE)
	{
		/* Variable size so call MKL with groups of size 1 */
		MKL_INT *group_size = (MKL_INT*) malloc(batch_count * sizeof(MKL_INT));
		for(int i = 0; i < batch_count; i++)
		{
			group_size[i] = 1;
		}
		/* Call MKL */
		cblas_zgemm_batch( BblasColMajor,
						   (CBLAS_TRANSPOSE*) transA, (CBLAS_TRANSPOSE*) transB,
						   (const MKL_INT*) M,
						   (const MKL_INT*) N,
						   (const MKL_INT*) K,
						   (const BBLAS_Complex64_t*) alpha,
						   (const BBLAS_MKL_Void_t**) arrayA, (const MKL_INT*) lda,
						   (const BBLAS_MKL_Void_t**) arrayB, (const MKL_INT*) ldb,
						   (const BBLAS_Complex64_t*) beta,
						   (BBLAS_MKL_Void_t**) arrayC, (const MKL_INT*) ldc,
						   batch_count, (const MKL_INT*) group_size);

	}else if (batch_opts == BBLAS_FIXED)
	{
		/* One large group */
		MKL_INT *group_size = (MKL_INT*) malloc(sizeof(MKL_INT));
		group_size[0] = batch_count;

		/* Call MKL */
		cblas_zgemm_batch( BblasColMajor,
						   (CBLAS_TRANSPOSE*) transA, (CBLAS_TRANSPOSE*) transB,
						   (const MKL_INT*) M,
						   (const MKL_INT*) N,
						   (const MKL_INT*) K,
						   (const BBLAS_Complex64_t*) alpha,
						   (const BBLAS_MKL_Void_t**) arrayA, (const MKL_INT*) lda,
						   (const BBLAS_MKL_Void_t**) arrayB, (const MKL_INT*) ldb,
						   (const BBLAS_Complex64_t*) beta,
						   (BBLAS_MKL_Void_t**) arrayC, (const MKL_INT*) ldc,
						   1, (const MKL_INT*) group_size);
	}else {
		/* batch_opts  is invalid */
		printf("ERROR: mkl_zgemm_batch called with invalid batch_opts parameter.\n");
	}
#endif
}
