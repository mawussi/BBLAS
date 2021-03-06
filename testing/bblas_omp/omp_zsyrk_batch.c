/**
 * @file omp_zsyrk_batch.c
 *
 * @brief BBLAS omp_zsyrk_batch double _Complex routine.
 *
 *  BBLAS is a software package provided by Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include<cblas.h>
#include "bblas_omp.h"
#include "bblas.h"
#include <omp.h>
#define COMPLEX

/**
    Purpose
    -------
    <b>zsyrk_batch</b> is an OpenMP version of zsyrk_batch.
    It performs the matrix-matrix operations

    arrayC[i] = alpha[i]*arrayA[i]*arrayA[i]**T + beta[i]*arrayC[i], or
	arrayC[i] = alpha[i]*arrayA[i]**T *arrayA[i] + beta[i]*arrayC[i],

    where alpha[i] and beta[i] are scalars, arrayC[i] is an N[i] by N[i] sym-
    metric matirix and arrayA[i] is N[i] by K[i] matrix in the
    first case and K[i] by N[i] matrix in the second case.


	Fixed and Variable Batch Operations
	-----------------------------------
	Two types of batch operation are supported depending upon the value of batch_opts.

	When <tt>batch_opts = BBLAS_VARIABLE</tt>
	- all parameters that are arrays must have length at least batch_count.
	- all parameters that are arrays must have all values set.

	When <tt>batch_opts = BBLAS_FIXED</tt>
	- all parameters that are arrays (except for arrayA, arrayC, and info)
	must have length at least one.
	- all parameters that are arrays (except for arrayA, arrayC, and info)
	need only to have their first value set.

	This means that for a <tt>BBLAS_FIXED</tt> batch,
	the values of uplo[0], trans[0], N[0], K[0],
	alpha[0], beta[0], lda[0], and ldc[0] are used for all computations.


    Parameters
    ----------
    @param[in]
    uplo    Array of <tt>enum BBLAS_UPLO</tt>.
            On entry, uplo[i] specifies whether the upper or
            lower triangular part of the matrix arrayC[i] is to
			be referenced as follows:
      -     = 'BblasUpper' Only the upper triangular part of
              the matrix is to be referenced.
      -     = 'BblasLower' Only the lower triangular part of
              the matrix is to be referenced.

    @param[in]
    trans   Array of <tt>enum BBLAS_TRANS</tt>.
            On entry, trans[i] specifies the operation to be
            performed as follows:
      -     = 'BblasNoTrans' arrayC[i] = alpha[i]*arrayA[i]*arrayA[i]**T + beta[i]*arrayC[i].
      -     = 'BblasTrans'   arrayC[i] = alpha[i]*arrayA[i]**T *arrayA[i] + beta[i]*arrayC[i].

    @param[in]
    N       Array of <tt>int</tt>.
            Each element N[i] specifies the number of rows and columns of the matrix
            arrayC[i]. N[i] must be greater than zero.

    @param[in]
    K       Array of <tt>int</tt>.
            On entry with trans[i] = 'BblasNoTrans', K[i] specifies the
            number of columns of the matrix arrayA[i],
			and upon entry with trans[i] = 'BblasTrans',
			K[i] specifies the number of rows of the matrix arrayA[i].
			K[i] must be greater than zero.

	@param[in]
    alpha   Array of <tt>complex_16</tt>.

    @param[in]
    arrayA   Array of pointers.
	         Each element arrayA[i] is a pointer to a COMPLEX_16 matrix of
			 dimension lda[i] by Ka[i],
			 where Ka[i] = K[i] when transA[i] = BblasNoTrans and is N[i] otherwise.
             Before entry with transA[i] = BblasNoTrans, the leading N[i] by K[i]
             part of the arrayA[i] must contain the elements of arrayA[i], otherwise
             the leading K[i] by N[i] part of the arrayA[i] must contain the
             elements of arrayA[i].

    @param[in]
    lda     Array of <tt>int</tt>.
            On entry, lda[i] specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When transA[i] = BblasNoTrans then
            lda[i] must be at least max( 1, N[i] ), otherwise lda[i] must be at
            least max( 1, K[i] ).

    @param[in]
    beta    Array of <tt>complex_16</tt>.
	        When beta[i] is set to zero arrayC[i] need not be set on input.

    @param[in,out]
    arrayC  Array of pointers.
	        Each elements arrayC[i] is a pointer to a COMPLEX_16 matrix of
			dimension ldc[i] by N[i].
            Before entry with uplo[i] = 'BblasUpper', the leading
            N[i] by N[i] upper triangular part of the arrayC[i] must con-
            tain the upper triangular part  of the symmetric
            matrix and the strictly lower triangular part of arrayC[i]
            is not referenced. On exit, the upper triangular
            part of the arrayC[i] is overwritten by the upper tri-
            angular part of the updated matrix.
			Before entry with uplo[i] = 'BlasLower', the leading N[i] by N[i] lower
            triangular part of the arrayC[i] must contain the lower
            triangular part  of the symmetric matrix  and the
            strictly upper triangular part of arrayC[i] is not refer-
            enced.  On exit, the lower triangular part of the
            arrayC[i] is overwritten by the lower triangular part
            of the updated matrix.

    @param[in]
    ldc     Array of <tt>int</tt>.
            On entry, ldc[i] specifies the first dimension of arrayC[i] as declared
            in the calling (sub) program. Each element ldc must be at least max( 1, N[i] ).

	@param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

    @param[in]
    batch_opts   <tt>enum BBLAS_OPTS</tt>
                 One of BBLAS_FIXED or BBLAS_VARIABLE depending upon the type of
				 batch operation required.

    @param[out]
    info   Array of <tt>int</tt>.
	       Each element info[i] is the error return code of the ith zsyrk in the batch,
		   these need not be set on entry.
  		   The error codes can be found in bblas_macros.h.
**/

void omp_zsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N,  const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, enum BBLAS_OPTS batch_opts, int *info)
{
/*Local variables  */
    int first_index = 0;
    int batch_iter;
    int LDA;
    char func_name[15] = "zsyrk_batch";

    /* Check input arguments */
    if (batch_count < 0)
    {
	xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
    }
    if (batch_opts == BBLAS_FIXED)
    {
        if ((uplo[first_index] != BblasUpper) &&
	    (uplo[first_index] != BblasLower))
	{
	    xerbla_batch(func_name, BBLAS_ERR_UPLO, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_UPLO;
            }
	    return;
	}
        if ((trans[first_index] != BblasNoTrans) &&
            (trans[first_index] != BblasTrans) &&
            (trans[first_index] != BblasConjTrans))
        {
            xerbla_batch(func_name, BBLAS_ERR_TRANS, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
                info[batch_iter]  = BBLAS_ERR_TRANS;
            }
            return;
        }
	if (N[first_index] < 0)
        {
	    xerbla_batch(func_name, BBLAS_ERR_N, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_N;
            }
	    return;
	}
        if (K[first_index] < 0)
        {
	    xerbla_batch(func_name, BBLAS_ERR_K, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_K;
            }
	    return;
	}
	if (trans[first_index] == BblasNoTrans)
        {
            LDA = N[first_index];
        } else
        {
            LDA = K[first_index];
        }
	if (lda[first_index] < max(1, LDA))
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDA, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] =  BBLAS_ERR_LDA;
            }
	    return;
	}
	if (ldc[first_index] < max(1, N[first_index]))
        {
	    xerbla_batch(func_name, BBLAS_ERR_LDC, first_index);
	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] = BBLAS_ERR_LDC;
            }
	    return;
	}
	/* particular case */
	if (N[first_index] == 0 || (K[first_index] == 0 ||
				    alpha[first_index] == (BBLAS_Complex64_t)0.0 ||
				    beta[first_index] == (BBLAS_Complex64_t)1.0))
	{
    	    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	    {
		info[batch_iter] =  BBLAS_SUCCESS;
            }
	    return;
	}
#pragma omp parallel for private(batch_iter)
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	{
            /*Call to cblas_zsyrk */
            cblas_zsyrk(
		BblasColMajor,
		uplo[first_index],
		trans[first_index],
		N[first_index],
		K[first_index],
		CBLAS_SADDR(alpha[first_index]),
		arrayA[batch_iter],
		lda[first_index],
		CBLAS_SADDR(beta[first_index]),
		arrayC[batch_iter],
		ldc[first_index]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
	} /*END FIXED SIZE FOR LOOP */
    }else if (batch_opts ==  BBLAS_VARIABLE)
    {

#pragma omp parallel for private(batch_iter, LDA)
        for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
        {
            /* Check input arguments */
            if ((uplo[batch_iter] != BblasUpper) &&
                (uplo[batch_iter] != BblasLower))
            {
                xerbla_batch(func_name, BBLAS_ERR_UPLO, batch_iter);
                info[batch_iter] = BBLAS_ERR_UPLO;
                continue;
            }
            if ((trans[batch_iter] != BblasNoTrans) &&
                (trans[batch_iter] != BblasTrans) &&
                (trans[batch_iter] != BblasConjTrans))
            {
                xerbla_batch(func_name, BBLAS_ERR_TRANS, batch_iter);
                info[batch_iter]  = BBLAS_ERR_TRANS;
                continue;
            }
            if (N[batch_iter] < 0)
            {
                xerbla_batch(func_name, BBLAS_ERR_N, batch_iter);
                info[batch_iter] = BBLAS_ERR_N;
                continue;
            }
            if (K[batch_iter] < 0)
            {
                xerbla_batch(func_name, BBLAS_ERR_K, batch_iter);
                info[batch_iter] = BBLAS_ERR_K;
                continue;
            }
            if (trans[batch_iter] == BblasNoTrans)
            {
                LDA = N[batch_iter];
            } else
            {
                LDA = K[batch_iter];
            }
            if (lda[batch_iter] < max(1, LDA))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDA, batch_iter);
                info[batch_iter] =  BBLAS_ERR_LDA;
                continue;
            }
            if (ldc[batch_iter] < max(1, N[batch_iter]))
            {
                xerbla_batch(func_name, BBLAS_ERR_LDC, batch_iter);
                info[batch_iter] = BBLAS_ERR_LDC;
                continue;
            }
            /* particular case */
            if (N[batch_iter] == 0  ||
		((K[batch_iter] == 0 || alpha[batch_iter] == (BBLAS_Complex64_t)0.0) &&
		 (beta[batch_iter] == (BBLAS_Complex64_t)1.0)))
            {
                info[batch_iter] =  BBLAS_SUCCESS;
                continue;
            }
            cblas_zsyrk(
		BblasColMajor,
		uplo[batch_iter],
		trans[batch_iter],
		N[batch_iter],
		K[batch_iter],
		CBLAS_SADDR(alpha[batch_iter]),
		arrayA[batch_iter],
		lda[batch_iter],
		CBLAS_SADDR(beta[batch_iter]),
		arrayC[batch_iter],
		ldc[batch_iter]);
            /* Successful */
            info[batch_iter] = BBLAS_SUCCESS;
        }
    }else
    {
	xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
    }
}
#undef COMPLEX
