/**
 * @file bblas_dutil.c
 *
 * @brief BBLAS testing utilities for double routines.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 * Contains routines used in the testing to modify randomly generated matrices
 * and compute the average error over an entire batch etc.
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from bblas_zutil.c normal z -> d, Mon Jun  6 09:44:14 2016
 **/
#endif

#include "bblas_common.h"
#if defined(BBLAS_WITH_MKL)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

/** Include real functions since using double real precision **/
#define REAL

/** Quick access to the matrix elements **/
#define A(i,j)  A[i + j*lda]

/**
 * Make a matrix symmetric/symmetric.  Makes diagonal real.
 * Sets Aji = conj( Aij ) for j < i, that is, copy & conjugate
 * lower triangle to upper triangle.
 **/

void bblas_dmake_symmetric(int lda, int N, double* A)
{
    int i, j;
    for( i=0; i < N; ++i ) {
        A(i,i) = creal( A(i,i) );
        for( j=0; j < i; ++j ) {
            A(j,i) = conj( A(i,j) );
        }
    }
}

#ifdef COMPLEX

/**
 * Make a matrix real-symmetric
 * Does NOT make diagonal real.
 * Sets Aji = Aij for j < i, that is,
 * copy lower triangle to upper triangle.
 **/

void bblas_dmake_symmetric(int lda, int N, double* A)
{
    int i, j;
    for( i=0; i < N; ++i ) {
        for( j=0; j < i; ++j ) {
            A(j,i) =  A(i,j);
        }
    }
}
#endif


/**
 * irandRange generates a random value (int)
 * in the range min_n and max_n
 **/

int irandRange(int min_n, int max_n)
{
    return rand() % (max_n - min_n + 1) + min_n;
}

/**
 * bblas_drandRange generates a random value (double)
 * in the range [0,max_n]: TODO replace by generic
 */

double bblas_drandRange( int max_n )
{
    return ( double )rand()/( double )( RAND_MAX/max_n );
}

/**
 * Computes statistics of the relative errors to summarise them for the user.
 **/

void bblas_dstatistic(bblas_dtest_t *test)
{

    enum BBLAS_ROUTINE routine  = test->routine;

    /*Compute avg(M), avg(N), avg(K) */
    if(test->batch_opts == BBLAS_VARIABLE)
    {
	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	    (routine == BBLAS_HEMM) || (routine == BBLAS_TRMM) ||
	    (routine == BBLAS_TRSM))
	{
	    test->avgM = bblas_avgarrayI(test->M, test->batch_count);
	}

	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYRK)  ||
	    (routine == BBLAS_HERK) || (routine == BBLAS_SYR2K) ||
	    (routine == BBLAS_HER2K))
	{
	    test->avgK = bblas_avgarrayI(test->K, test->batch_count);
	}

	test->avgN = bblas_avgarrayI(test->N, test->batch_count);

    } else if (test->batch_opts == BBLAS_FIXED)
    {

	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	    (routine == BBLAS_HEMM) || (routine == BBLAS_TRMM) ||
	    (routine == BBLAS_TRSM))
	{
	    test->avgM = test->M[0];
	}

	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYRK)  ||
	    (routine == BBLAS_HERK) || (routine == BBLAS_SYR2K) ||
	    (routine == BBLAS_HER2K))
	{
	    test->avgK = test->K[0];
	}
	test->avgN = test->N[0];

    } else
    {
	bblas_error("testing_dgemm_batch.c", "wrong batch_opts value");
    }

    /*Statistics on the error */
    switch(test->target)
    {
    case BBLAS_MKL:
      test->mkl_min_error = bblas_dminarrayD(test->mkl_error, test->batch_count);
      test->mkl_avg_error = bblas_davgarrayD(test->mkl_error, test->batch_count);
      test->mkl_max_error = bblas_dmaxarrayD(test->mkl_error, test->batch_count);
      test->mkl_std_error = bblas_dstdarrayD(test->mkl_error, test->batch_count);
      break;

    case BBLAS_CUBLAS:
    case BBLAS_MAGMA:
      test->device_min_error = bblas_dminarrayD(test->device_error, test->batch_count);
      test->device_avg_error = bblas_davgarrayD(test->device_error, test->batch_count);
      test->device_max_error = bblas_dmaxarrayD(test->device_error, test->batch_count);
      test->device_std_error = bblas_dstdarrayD(test->device_error, test->batch_count);
      break;

    case BBLAS_OTHER:
      test->other_min_error = bblas_dminarrayD(test->other_error, test->batch_count);
      test->other_avg_error = bblas_davgarrayD(test->other_error, test->batch_count);
      test->other_max_error = bblas_dmaxarrayD(test->other_error, test->batch_count);
      test->other_std_error = bblas_dstdarrayD(test->other_error, test->batch_count);
      break;

    default:
      printf("In bblas_dstatistic(): Target no defined\n");
      exit(EXIT_FAILURE);
    }
}


/**
 * Print a matrix.
 **/
void bblas_dprintmatrix(double *matrix, int row, int col)
{
    /*Local variables */
    int i, j;

    for (i=0; i < row; i++)
    {
	printf("\n\n");
	for (j=0; j < col; j++)
	{
#ifdef COMPLEX
	    printf("%1.2f + %1.2f\t", creal(matrix[i*col+j]), cimag(matrix[i*col+j]));
#else
	    printf("%1.2f",matrix[i*col+j]);
#endif
	}
    }
    printf("\n");
}

/**
 * Decide whether a batch is fixed or variable.
 **/

char* bblas_getoption(enum BBLAS_OPTS opts)
{
   /*Local variable */
    char funcname[] = "bblas_getoption";

    switch(opts)
    {
	case BBLAS_VARIABLE:
	    return "BATCH OPTION: VARIABLE";
	    break;

	case BBLAS_FIXED:
	    return "BATCH OPTION: FIXED";
	    break;

	default:
	    printf("ERROR in %s, undefined bblas routine name\n",funcname);
	    exit(EXIT_FAILURE);
    }

}

/**
 * Get the name of the current test routine.
 **/

char* bblas_getroutine(enum BBLAS_ROUTINE routine)
{
    /*Local variable */
    char funcname[] = "bblas_getroutine";

    switch(routine)
    {
	case BBLAS_GEMM:
	    return "DGEMM";
	    break;

	case BBLAS_HEMM:
	    return "DSYMM";
	    break;

	case BBLAS_HER2K:
	    return "DSYR2K";
	    break;

	case BBLAS_HERK:
	    return "DSYRK";
	    break;

	case BBLAS_SYMM:
	    return "DSYMM";
	    break;

	case BBLAS_SYR2K:
	    return "DSYR2K";
	    break;

	case BBLAS_SYRK:
	    return "DSYRK";
	    break;

	case BBLAS_TRMM:
	    return "DTRMM";
	    break;

	case BBLAS_TRSM:
	    return "DTRSM";
	    break;

	default:
	    printf("ERROR in %s, undefined bblas routine name\n",funcname);
	    exit(EXIT_FAILURE);
    }
}



/**
 * Computes the maximum value of an array of double precision floats.
 **/

double bblas_dmaxarrayD(double *myArray, int size)
{
    int iter;
    double maxValue = myArray[0];

    for (iter = 0; iter < size; ++iter)
    {
        if ( myArray[iter] > maxValue )
	{
            maxValue = myArray[iter];
        }
    }
    return maxValue;
}


/**
 * Computes the minimum value of an array of double precision floats.
 **/

double bblas_dminarrayD(double *myArray, int size)
{
    int iter;
    double minValue = myArray[0];

    for (iter = 0; iter < size; ++iter)
    {
        if ( myArray[iter] < minValue )
	{
            minValue = myArray[iter];
        }
    }
    return minValue;
}

/**
 * Computes the mean value of an array of double precision floats.
 **/

double bblas_davgarrayD(double *myArray, int size)
{
    int iter;
    double avg = 0.;

    for (iter = 0; iter < size; ++iter)
    {
	avg += myArray[iter];
    }
    return avg/size;
}

/**
 * Computes the standard deviation of an array of double precision floats.
 **/

double bblas_dstdarrayD(double *myArray, int size)
{
    int iter;
    double avg, sd=0.;

    avg =  bblas_davgarrayD(myArray, size);

    for (iter = 0; iter < size; ++iter)
    {
	sd += (myArray[iter] -avg)*(myArray[iter] -avg);
    }
    return sd/size;
}

/**
 * Computes the minimum value of an array of integers.
 **/

int bblas_minarrayI(int *myArray, int size)
{
    int iter;
    int minValue = myArray[0];

    for (iter = 0; iter < size; ++iter)
    {
        if ( myArray[iter] < minValue )
	{
            minValue = myArray[iter];
        }
    }
    return minValue;
}

/**
 * Computes the mean of an array of integers.
 **/

int bblas_avgarrayI(int *myArray, int size)
{
    int iter;
    int avg = 0;

    for (iter = 0; iter < size; ++iter)
    {
	avg += myArray[iter];
    }
    return avg/size;
}


/**
 * Transform BBLAS enum values for <tt>trans</tt>, <tt>uplo</tt> etc. to human-readable strings.
 **/

char* bblas_op2char(unsigned int op)
{
    char *opname = (char*)malloc(30*sizeof(char));
    switch(op)
    {
        case BblasNoTrans:
            strcpy(opname,"CblasNoTrans");
            break;

        case BblasTrans:
	    strcpy(opname,"CblasTrans");
            break;

        case BblasConjTrans:
            strcpy(opname,"CblasConjTrans");
            break;

        case BblasLower:
            strcpy(opname,"CblasLower");
            break;

        case BblasUpper:
	    strcpy(opname,"CblasUpper");
            break;

        case BblasNonUnit:
	    strcpy(opname,"CblasNonUnit");
            break;

        case BblasUnit:
	   strcpy(opname,"CblasUnit");
            break;

        case BblasLeft:
	   strcpy(opname,"CblasLeft");
            break;

        case BblasRight:
	    strcpy(opname,"CblasRight");
            break;

        default:
            return 0;
            exit(EXIT_FAILURE);
    }
    return opname;
}

/**
 * Get the amount of data needed.
 **/

int bblas_dnbdata(bblas_dtest_t *test)
{
    enum BBLAS_OPTS batch_opts        = test->batch_opts;
    int nb_data;
    char function_name[NAME_LENGTH]   ="bblas_getdatacount";

    switch(batch_opts)
    {
	case BBLAS_VARIABLE:
	    nb_data =  test->batch_count;
	    break;

	case BBLAS_FIXED:
	    nb_data = 1;
	    break;

	default:
	    bblas_fatal_error(function_name, "wrong batch_opts value");
    }

    return nb_data;
}


/**
 * Inject an error into the computation if this is set in the input file.
 **/

void bblas_dset_error(bblas_dtest_t *test)
{
    int   nb_data =  bblas_dnbdata(test);
    int   routine = test->routine;
    int   error_index = irandRange(0, nb_data);

    if (test->global_error){
	test->batch_count = -1;
	return;
    }


    if (test->batch_opts == BBLAS_FIXED)
    {
	error_index = 0;
    }

    if (test->set_error)
    {
	if ((routine == BBLAS_GEMM) || (routine == BBLAS_SYMM) ||
	    (routine == BBLAS_HEMM) || (routine == BBLAS_TRMM) ||
	    (routine == BBLAS_TRSM))
	{
	    test->M[error_index] = -1;
	}else if ((routine == BBLAS_SYRK) || (routine == BBLAS_HERK) ||
		  (routine == BBLAS_SYR2K)|| (routine == BBLAS_HER2K))
	{
	    test->K[error_index] = -1;
	}
    }
}

/**
 * Check whether a computation has passed or failed our accuracy test.
 * When new_accuracy=1 in the input file this uses an appropriate
 * forward/backward error bound,
 * otherwise this looks at the relative error.
 **/

void bblas_dpassed_failed(bblas_dtest_t  *test, double error, char *result, int info)
{
    double eps = LAPACKE_dlamch_work('e');

	/* Use our new accuracy test based on forward/backward error analysis*/
    if (test->new_accuracy)
    {
	if( (error > 1) || (info))
	{
	    strcpy(result, "FAILED");

	}else
	{
	    strcpy(result, "PASSED");
	}
    }else
    {
		/* Use old accuracy test based on the relative error */
	if((error > eps*test->tolerance) || (info))
	{
	    strcpy(result, "FAILED");

	}else
	{
	    strcpy(result, "PASSED");
	}
    }
}

/**
 * Set the batch_count.
 **/

void bblas_dset_batch_count(bblas_dtest_t *test)
{

  test->batch_count = test->minbatch_count* (test->current_iter+1);
}


#undef REAL
