/**
 * @file bblas_saccuracy.c
 *
 * @brief BBLAS accuracy tests for float_Complex precision.
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
 * @generated from bblas_zaccuracy.c normal z -> s, Mon Jun  6 09:44:14 2016
 **/
#endif

#include "bblas_common.h"
#include <lapacke.h>
#include <cblas.h>



#define REAL

/**
 * bblas_sgemm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_{\infty} +
 *           |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_sgemm_tolerance(bblas_stest_t *test, float **C_final)
{

    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, Bm, Bn, M, K, N, ldc;
    int max_work_size = max(test->maxK,max(test->maxM, test->maxN));
    float Anorm, Bnorm, Rnorm, result;
    float eps = LAPACKE_slamch_work('e');
    float alpha =-1;
    float *work = (float *)malloc(max_work_size*sizeof(float));

    /*Temporary buffer to save (test-arrayC -C_final) */
    float **C_diff;
    C_diff = (float **) malloc(batch_count*sizeof(float *));

    /*Make a copy of test->arrayC in C_diff */
    bblas_scopy_Cinit(test, C_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    M = test->M[batch_iter];
	    K   = test->K[batch_iter];
	    N   = test->N[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->transA[batch_iter] == BblasNoTrans) {
		Am = M; An = K;
	    } else {
		Am = K; An = M;
	    }
	    if (test->transB[batch_iter] == BblasNoTrans) {
		Bm = K; Bn = N;
	    } else {
		Bm = N; Bn = K;
	    }

	    /*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', M, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    result = Rnorm / ((K*cabs(test->alpha[batch_iter])*Anorm*Bnorm + cabs(test->beta[batch_iter])*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	M   = test->M[first_index];
	K   = test->K[first_index];
	ldc = test->ldc[first_index];
	N   = test->N[first_index];

	if (test->transA[first_index] == BblasNoTrans) {
	    Am = M; An = K;
	} else {
	    Am = K; An = M;
	}
	if (test->transB[first_index] == BblasNoTrans) {
	    Bm = K; Bn = N;
	} else {
	    Bm = N; Bn = K;
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
		/*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);
	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', M, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the relative error */
	    result = Rnorm / ((K*cabs(test->alpha[first_index])*Anorm*Bnorm + cabs(test->beta[first_index])*test->Cinitnorm[batch_iter])*eps);

	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal: Fixed, Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_sstarmm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_{\infty} +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_sstarmm_tolerance(bblas_stest_t *test, float **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, Bm, Bn, N, ldc;
    int max_work_size = max(test->maxM, test->maxN);
    float Anorm, Bnorm, Rnorm, result;
    float eps = LAPACKE_slamch_work('e');
    float alpha =-1;
    float *work = (float *)malloc(max_work_size*sizeof(float));

    /*Temporary buffer to save (test-arrayC -C_final) */
    float **C_diff;
    C_diff = (float **) malloc(batch_count*sizeof(float *));

    /*Make a copy of test->arrayC in C_diff */
    bblas_scopy_Cinit(test, C_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		Am = test->M[batch_iter];
	    } else {
		Am = test->N[batch_iter];
	    }
	    Bm = test->ldb[batch_iter];
	    Bn = test->N[batch_iter];
	    /*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, Am, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    result = Rnorm / ((N*cabs(test->alpha[batch_iter])*Anorm*Bnorm + cabs(test->beta[batch_iter])*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	N = test->N[first_index];
	ldc = test->ldc[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    Am = test->M[first_index];
	} else {
	    Am = test->N[first_index];
	}
	Bm = test->ldb[first_index];
	Bn = test->N[first_index];

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, Am, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the relative error */
	    result = Rnorm / ((N*cabs(test->alpha[first_index])*Anorm*Bnorm + cabs(test->beta[first_index])*test->Cinitnorm[batch_iter])*eps);

	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal: Fixed, Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_sstark_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|A\|_1 +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_sstark_tolerance(bblas_stest_t *test, float **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, N, K, ldc;
    char u;
    int max_work_size = max(test->maxK, test->maxN);
    float Anorm, ATnorm, Rnorm, result;
    float eps = LAPACKE_slamch_work('e');
    float alpha_norm, beta_norm;
    float alpha =-1;
    float *work = (float *)malloc(max_work_size*sizeof(float));

    /*Temporary buffer to save (test-arrayC -C_final) */
    float **C_diff;
    C_diff = (float **) malloc(batch_count*sizeof(float *));


    /*Make a copy of test->arrayC in C_diff */
    bblas_scopy_Cinit(test, C_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    K = test->K[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->trans[batch_iter] == BblasNoTrans) {
		Am = N; An = K;
	    } else {
		Am = K; An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	      {
		u = 'L';
	      }else
	      {
		u = 'U';
	      }
	    /*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

	    /*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slansy_work(LAPACK_COL_MAJOR, 'I', u, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the one norm of A */
	    ATnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

	    if (test->routine == BBLAS_HERK)
	    {
		beta_norm   = cabs(test->beta_herk[batch_iter]);
		alpha_norm  = cabs(test->alpha_herk[batch_iter]);
	    }else{
		beta_norm  = cabs(test->beta[batch_iter]);
		alpha_norm  = cabs(test->alpha[batch_iter]);
	    }

	    result = Rnorm / ((K*alpha_norm*Anorm*ATnorm + beta_norm*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	    N = test->N[first_index];
	    K = test->K[first_index];
	    ldc = test->ldc[first_index];

	    if (test->trans[first_index] == BblasNoTrans) {
		Am = N; An = K;
	    } else {
		Am = K; An = N;
	    }
	    if(test->uplo[first_index] == BblasLower )
	      {
		u = 'L';
	      }else
	      {
		u = 'U';
	      }

	    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	      {
		/*Compute the error C - C_finial */
		cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		Rnorm = LAPACKE_slansy_work(LAPACK_COL_MAJOR, 'I', u, N, C_diff[batch_iter], ldc, work);

		/*Compute the infinity norm of A */
		Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the one norm of A */
		ATnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

		if (test->routine == BBLAS_HERK)
		{
		    beta_norm   = cabs(test->beta_herk[first_index]);
		    alpha_norm  = cabs(test->alpha_herk[first_index]);
		}else{
		    beta_norm   = cabs(test->beta[first_index]);
		    alpha_norm  = cabs(test->alpha[first_index]);
		}

		/*Compute the relative error */
		result = Rnorm / ((K*alpha_norm*Anorm*ATnorm +
				   beta_norm*test->Cinitnorm[batch_iter])*eps);

		switch(test->target)
		{
		    case BBLAS_MKL:
			test->mkl_error[batch_iter] = result;
			break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		  test->device_error[batch_iter] = result;
		  break;

		    case BBLAS_OTHER:
			test->other_error[batch_iter] = result;
			break;
		    default:
			printf("In bblas_scheck_Cfinal: Fixed, Target no defined\n");
			exit(EXIT_FAILURE);
		}
	    }
    }else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_sstar2k_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((K|\alpha|\|A\|_{\infty}\|B\|_1 +
 * K|\alpha|\|B\|_{\infty}\|A\|_1  +
 * |\beta|\|C_{\mathrm{init}}\|_{\infty})\epsilon).@f]
 *
 **/

void bblas_sstar2k_tolerance(bblas_stest_t *test, float **C_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, Bn, Bm, N, K, ldc;
    int max_work_size = max(test->maxK, test->maxN);
    float Anorm, Bnorm, ATnorm, BTnorm,  Rnorm, result;
    float beta_norm;
    float eps = LAPACKE_slamch_work('e');
    float alpha =-1;
    float *work = (float *)malloc(max_work_size*sizeof(float));

    /*Temporary buffer to save (test-arrayC -C_final) */
    float **C_diff;
    C_diff = (float **) malloc(batch_count*sizeof(float *));

    /*Make a copy of test->arrayC in C_diff */
    bblas_scopy_Cinit(test, C_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    K = test->K[batch_iter];
	    ldc = test->ldc[batch_iter];

	    if (test->trans[batch_iter] == BblasNoTrans) {
		Am = N; An = K;
		Bn = K;
	    } else {
		Am = K; An = N;
		Bn = N;
	    }
	    Bm = test->ldb[batch_iter];

	    /*Compute the error C - C_finial */
	    cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

	    /*Compute the infinity norm of A */
	    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the one norm of A */
	    ATnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

	    /*Compute the infinity norm of B */
	    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    /*Compute the one norm of B */
	    BTnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Bm, Bn, test->arrayB[batch_iter], Bm, work);

	    if(test->routine == BBLAS_HER2K)
	    {
		beta_norm = cabs(test->beta_herk[batch_iter]);
	    }else
	    {
		beta_norm = cabs(test->beta[batch_iter]);
	    }

	    result = Rnorm / ((K*cabs(test->alpha[batch_iter])*Anorm*BTnorm +
			       K*cabs(test->alpha[batch_iter])*Bnorm*ATnorm +
			       beta_norm*test->Cinitnorm[batch_iter])*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}

    }else  if( test->batch_opts == BBLAS_FIXED ) //fixed  size
    {
	N = test->N[first_index];
	K = test->K[first_index];
	ldc = test->ldc[first_index];

	if (test->trans[first_index] == BblasNoTrans) {
		Am = N; An = K;
		Bn = K;
	} else {
	    Am = K; An = N;
	    Bn = N;
	}
	Bm = test->ldb[first_index];


	    for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	    {
		/*Compute the error C - C_finial */
		cblas_saxpy (ldc*N, (alpha), C_final[batch_iter], 1, C_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
		Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldc, N, C_diff[batch_iter], ldc, work);

		/*Compute the infinity norm of A */
		Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the one norm of A */
		ATnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Am, An, test->arrayA[batch_iter], Am, work);

		/*Compute the infinity norm of B */
		Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', Bm, Bn, test->arrayB[batch_iter], Bm, work);

		/*Compute the one norm of B */
		BTnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'O', Bm, Bn, test->arrayB[batch_iter], Bm, work);

		if(test->routine == BBLAS_HER2K)
		{
		    beta_norm = cabs(test->beta_herk[first_index]);
		}else
		{
		    beta_norm = cabs(test->beta[first_index]);
		}

		/*Compute the relative error */
		result = Rnorm / ((K*cabs(test->alpha[first_index])*Anorm*BTnorm +
				   K*cabs(test->alpha[first_index])*Bnorm*ATnorm +
				   beta_norm*test->Cinitnorm[batch_iter])*eps);


		switch(test->target)
		{
		case BBLAS_MKL:
		  test->mkl_error[batch_iter] = result;
		  break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		  test->device_error[batch_iter] = result;
		  break;

		case BBLAS_OTHER:
		  test->other_error[batch_iter] = result;
		  break;
		default:
		  printf("In bblas_scheck_Cfinal: Fixed, Target no defined\n");
		  exit(EXIT_FAILURE);
		}
	    }
    }else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

/*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(C_diff[index]);
    }
    free(C_diff);
    free(work);
}


/**
 * bblas_strmm_tolerance computes the forward error
 * @f[\|C - C_{\mathrm{seq}}\|_{\infty}/((|\alpha|\|A\|_{\infty}\|B\|_{\infty}N\epsilon).@f]
 *
 **/

void bblas_strmm_tolerance(bblas_stest_t *test, float **B_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, ldb, N, M;
    int max_work_size = max(test->maxK, test->maxN);
    char u, diag;
    float Anorm, Rnorm, result;
    float eps = LAPACKE_slamch_work('e');
    float alpha =-1;
    float *work = (float *)malloc(max_work_size*sizeof(float));



    /*Temporary buffer to save (test-arrayB -B_final) */
    float **B_diff;
    B_diff = (float **) malloc(batch_count*sizeof(float *));

    /*Make a copy of test->arrayB in B_diff */
    bblas_scopy_Binit(test, B_diff);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    M = test->M[batch_iter];
	    ldb = test->ldb[batch_iter];

	    Am = test->lda[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		An = M;
	    }else
	    {
		An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	    {
	       u = 'L';
	    }else
	    {
	       u = 'U';
	    }

	    if(test->diag[batch_iter] == BblasUnit )
	    {
		diag = 'U';
	    }else
	    {
		diag = 'N';
	    }

	    /*Compute the error B - B_finial */
	    cblas_saxpy (ldb*N, (alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, B_diff[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_slantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);
	    
		result = Rnorm/(N*cabs(test->alpha[batch_iter])*Anorm*test->Binitnorm[batch_iter]*eps);
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }

	}
    }else  if( test->batch_opts == BBLAS_FIXED )
    {

	N = test->N[first_index];
	M = test->M[first_index];
	ldb = test->ldb[first_index];

	Am = test->lda[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    An = M;
	}else
	{
	    An = N;
	}

	if(test->uplo[first_index] == BblasLower )
	{
	    u = 'L';
	}else
	{
	    u = 'U';
	}

	if(test->diag[first_index] == BblasUnit )
	{
	    diag = 'U';
	}else
	{
	    diag = 'N';
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    /*Compute the error B - B_finial */
	    cblas_saxpy (ldb*N, (alpha), B_final[batch_iter], 1, B_diff[batch_iter], 1);

		/*Compute  the infinity norm associated  with the error */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, B_diff[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_slantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    result = Rnorm/(N*cabs(test->alpha[first_index])*Anorm*test->Binitnorm[batch_iter]*eps);

	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }
    else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

    /*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(B_diff[index]);
    }
    free(B_diff);
    free(work);
}

/**
 * bblas_strsm_tolerance computes the backward error
 * @f[\|\alpha B - AX \|_{\infty}/((\|A\|_{\infty}\|X\| + |\alpha|\|B\|_{\infty})N\epsilon).@f]
 *
 **/

void bblas_strsm_tolerance(bblas_stest_t *test, float **B_final)
{
    /*Local variables */
    int batch_count = test->batch_count;
    int batch_iter, first_index = 0;
    int Am, An, ldb, N, M;
    int max_work_size = max(test->maxK, test->maxN);
    char u, diag;
    float Anorm, Rnorm, result, Xnorm;
    float eps = LAPACKE_slamch_work('e');
    float alpha;
    float *work = (float *)malloc(max_work_size*sizeof(float));



    /*Variable to save the residual alphaB -AX */
    float **residual;
    residual = (float **) malloc(batch_count*sizeof(float *));

    /*Make a copy of test->arrayB (X) in residual */
    bblas_scopy_Binit(test, residual);

    /*Variable size */
    if( test->batch_opts == BBLAS_VARIABLE )
    {
	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{
	    N = test->N[batch_iter];
	    M = test->M[batch_iter];
	    ldb = test->ldb[batch_iter];

	    Am = test->lda[batch_iter];

	    if(test->side[batch_iter] == BblasLeft )
	    {
		An = M;
	    }else
	    {
		An = N;
	    }

	    if(test->uplo[batch_iter] == BblasLower )
	    {
	       u = 'L';
	    }else
	    {
	       u = 'U';
	    }

	    if(test->diag[batch_iter] == BblasUnit )
	    {
		diag = 'U';
	    }else
	    {
		diag = 'N';
	    }

	    /*Compute residual <- AX  (strmm) */
	    alpha = 1;
	    cblas_strmm(BblasColMajor, test->side[batch_iter], test->uplo[batch_iter], test->transA[batch_iter],
			test->diag[batch_iter], M, N, (alpha), test->arrayA[batch_iter],
			test->lda[batch_iter], residual[batch_iter], ldb);

	    /*Compute Binit = alpha Binit */
	    for (int index=0;  index < ldb*N; index++)
	    {
		test->Binit[batch_iter][index] = test->alpha[batch_iter]*test->Binit[batch_iter][index];
	    }

	    /*Compute the  residual <- alpha B - Ax  */
	    alpha= -1;
	    cblas_saxpy (ldb*N, (alpha), test->Binit[batch_iter], 1, residual[batch_iter], 1);
		/*Compute  the infinity norm associated  with the residual */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, residual[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_slantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    /*Compute  the infinity norm associated  with X (B) */
	    Xnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, test->arrayB[batch_iter], ldb, work);

	    result = Rnorm/(Anorm*Xnorm*N*eps);
	    
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }

	}
    }else  if( test->batch_opts == BBLAS_FIXED )
    {

	N = test->N[first_index];
	M = test->M[first_index];
	ldb = test->ldb[first_index];

	Am = test->lda[first_index];

	if(test->side[first_index] == BblasLeft )
	{
	    An = M;
	}else
	{
	    An = N;
	}

	if(test->uplo[first_index] == BblasLower )
	{
	    u = 'L';
	}else
	{
	    u = 'U';
	}

	if(test->diag[first_index] == BblasUnit )
	{
	    diag = 'U';
	}else
	{
	    diag = 'N';
	}

	for( batch_iter =0; batch_iter < batch_count ; batch_iter++)
	{

	    /*Compute residual <- AX  (strmm) */
	    alpha = 1;
	    cblas_strmm(BblasColMajor, test->side[first_index], test->uplo[first_index], test->transA[first_index],
			test->diag[first_index], M, N, (alpha), (const float*) test->arrayA[batch_iter],
			test->lda[first_index], residual[batch_iter], ldb);

	    /*Compute Binit = alpha Binit */
	    for (int index=0;  index < ldb*N; index++)
	    {
		test->Binit[batch_iter][index] = test->alpha[first_index]*test->Binit[batch_iter][index];
	    }

	    /*Compute the  residual <- alpha B - Ax  */
	    alpha = -1;
	    cblas_saxpy (ldb*N, (alpha), test->Binit[batch_iter], 1, residual[batch_iter], 1);

		/*Compute  the infinity norm associated  with the residual */
	    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, residual[batch_iter], ldb, work);

	    /*Infinity norm of the triangular matrix A */
	    Anorm = LAPACKE_slantr_work(LAPACK_COL_MAJOR, 'I', u, diag, Am,
					An, test->arrayA[batch_iter], Am, work);

	    /*Compute  the infinity norm associated  with X (B) */
	    Xnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'I', ldb, N, test->arrayB[batch_iter], ldb, work);

	    result = Rnorm/(Anorm*Xnorm*N*eps);
	    
	    /*Compute the relative error */
	    switch(test->target)
	    {
		case BBLAS_MKL:
		    test->mkl_error[batch_iter] = result;
		    break;

		case BBLAS_CUBLAS:
		case BBLAS_MAGMA:
		    test->device_error[batch_iter] = result;
		    break;

		case BBLAS_OTHER:
		    test->other_error[batch_iter] = result;
		    break;
		default:
		    printf("In bblas_scheck_Cfinal, Variable: Target no defined\n");
		    exit(EXIT_FAILURE);
	    }
	}
    }
    else
    {
	bblas_error("bblas_stesting.c", "wrong batch_opts value");
    }

    /*Free allocated memory */
    for(int index =0; index < batch_count; index++)
    {
	free(residual[index]);
    }
    free(residual);
    free(work);
}


/**
 * Calls bblas_sstar2k_tolerance.
 **/
#ifdef COMPLEX
void bblas_ssyr2k_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstar2k_tolerance(test, C_final);

}
#endif
/**
 * Calls bblas_sstar2k_tolerance.
 **/

void bblas_ssyr2k_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstar2k_tolerance(test, C_final);
}

/**
 * Calls bblas_sstark_tolerance.
 **/

#ifdef COMPLEX
void bblas_ssyrk_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstark_tolerance(test, C_final);
}
#endif
/**
 * Calls bblas_sstark_tolerance.
 **/

void bblas_ssyrk_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstark_tolerance(test, C_final);
}

/**
 * Calls bblas_sstarmm_tolerance.
 **/

#ifdef COMPLEX
void bblas_ssymm_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstarmm_tolerance(test, C_final);
}
#endif
/**
 * Calls bblas_sstamm_tolerance.
 **/

void bblas_ssymm_tolerance(bblas_stest_t *test, float **C_final)
{
    bblas_sstarmm_tolerance(test, C_final);
}
