/**
 * @file testing_strmm_batch.c
 *
 * @brief Testing strmm_batch float routine.
 *
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-02-20
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from testing_ztrmm_batch.c normal z -> s, Mon Jun  6 09:44:14 2016
 **/
#endif

#define REAL

#include "bblas_common.h"
#include "bblas_mkl.h"
#include "bblas.h"
#include "bblas_omp.h"

/**
 * Makes calls to the reference implementation, CuBLAS, MKL, etc.
 * gets the statistics of the associated errors and prints the results.
 **/

void testing_strmm_batch(bblas_stest_t *test){

    for(test->current_iter = 0; test->current_iter < test->nb_test; ++test->current_iter )
    {
	/* set test parameters */
	bblas_ssettest(test);

	/* Copy the arrayB in temporay buffer for accuracy testing */
	switch (test->target)
	{
	    case BBLAS_MKL:
		bblas_scopy_Binit(test, test->mkl_result);
		break;

	    case BBLAS_CUBLAS:
		bblas_scopy_Binit(test, test->device_result);
		break;

	    case BBLAS_OTHER:
		bblas_scopy_Binit(test, test->other_result);
		break;

	    case BBLAS_CUMKL:
		bblas_scopy_Binit(test, test->mkl_result);
		bblas_scopy_Binit(test, test->device_result);
		break;

	    default:
		printf("In testing_strmm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Compute flops */
	bblas_sflops(test);

	/* Error injection */
	if (test->current_iter == test->faulty_iter)
	{
	    bblas_sset_error(test);
	}


	/* =====================================================================
	   Performs operation using SEC
	   =================================================================== */
	test->ref_time = bblas_wtime();

	bblas_strmm_batch(test->side, test->uplo, test->transA, test->diag,
		    test->M, test->N, test->alpha, (const float**)test->arrayA,
		    test->lda, test->arrayB, test->ldb,	test->batch_count,
		    test->batch_opts, test->info);

	test->ref_time = bblas_wtime() - test->ref_time;

	/*Compute seq_perf */
	test->ref_perf = (test->flops/1e9)/test->ref_time;

	switch(test->target)
	{
	    case BBLAS_MKL:
		/* =====================================================================
		   Performs operation using  MKL BBLAS
		   =================================================================== */

		test->mkl_time = bblas_wtime();
		bblas_strmm_batch(test->side, test->uplo, test->transA, test->diag,
			    test->M, test->N, test->alpha, (const float**)test->arrayA,
			    test->lda, test->mkl_result, test->ldb, test->batch_count,
			    test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_strmm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_scheck_Bfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*Print results*/
		bblas_sprint_result(test);

		break;

	    case BBLAS_CUBLAS:
		/* =====================================================================
		   Performs operation using CUBLAS BBLAS
		   =================================================================== */

		test->device_time = bblas_wtime();
		bblas_strmm_batch(test->side, test->uplo, test->transA, test->diag,
			    test->M, test->N, test->alpha, (const float**)test->arrayA,
			    test->lda, test->device_result, test->ldb,	test->batch_count,
			    test->batch_opts, test->info);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_strmm_tolerance(test, test->device_result);
		} else
		{
		    bblas_scheck_Bfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*Print results*/
		bblas_sprint_result(test);

		break;

	    case BBLAS_CUMKL:

		/*******************************
		 *        [1] CUBLAS CASE      *
		********************************/

		/*Change target to particular case of BBLAS_CUBLAS */
		test->target = BBLAS_CUBLAS;

		test->device_time = bblas_wtime();
		bblas_strmm_batch(test->side, test->uplo, test->transA, test->diag,
			    test->M, test->N, test->alpha, (const float**)test->arrayA,
			    test->lda, test->device_result, test->ldb, test->batch_count,
			    test->batch_opts, test->info);

		test->device_time = bblas_wtime() - test->device_time;

		/*Compute perf */
		test->device_perf = (test->flops/1e9)/test->device_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_strmm_tolerance(test, test->device_result);
		} else
		{
		    bblas_scheck_Bfinal(test, test->device_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*******************************
		 *        [2] MKL CASE      *
		********************************/

		/*Change target to particular case of BBLAS_MKL */
		test->target = BBLAS_MKL;

		test->mkl_time = bblas_wtime();
		bblas_strmm_batch(test->side, test->uplo, test->transA, test->diag,
			    test->M, test->N, test->alpha, (const float**)test->arrayA,
			    test->lda, test->mkl_result, test->ldb, test->batch_count,
			    test->batch_opts, test->info);

		test->mkl_time = bblas_wtime() - test->mkl_time;

		/*Compute perf */
		test->mkl_perf = (test->flops/1e9)/test->mkl_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_strmm_tolerance(test, test->mkl_result);
		} else
		{
		    bblas_scheck_Bfinal(test, test->mkl_result);
		}

		/*Perform a statistic */
		bblas_sstatistic(test);

		/*******************************
		 *[3] BOTH CUBLAS & MKL CASE   *
		********************************/

		/*Change target to BBLAS_CUMKL */
		test->target = BBLAS_CUMKL;

		/*Print results*/
		bblas_sprint_result(test);
		break;

	    case BBLAS_OTHER:
		/* =====================================================================
		   Performs operation using OMP CPU BLAS
		   =================================================================== */

		test->other_time = bblas_wtime();
		omp_strmm_batch(test->side, test->uplo, test->transA, test->diag,
			    test->M, test->N, test->alpha, (const float**)test->arrayA,
			    test->lda, test->other_result, test->ldb,	test->batch_count,
			    test->batch_opts, test->info);

		test->other_time = bblas_wtime() - test->other_time;


		/*Compute perf */
		test->other_perf = (test->flops/1e9)/test->other_time;

		/*Compute relative errors */
		if(test->new_accuracy)
		{
		    bblas_strmm_tolerance(test, test->other_result);
		} else
		{
		    bblas_scheck_Bfinal(test, test->other_result);
		}
		/*Perform a statistic */
		bblas_sstatistic(test);

		/*Print results*/
		bblas_sprint_result(test);
		break;

	    default:
		printf("In testing_strmm_batch(): Target no defined\n");
		exit(EXIT_FAILURE);
	}

	/*Free batch configuration */

	bblas_sfreetest(test);
	printf("\n");
    }
}


#undef REAL
