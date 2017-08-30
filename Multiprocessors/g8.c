/*****************************************************
 *
 * Gaussian elimination
 *
 * parallel version
 *
 *****************************************************/

#include <stdio.h>
#include <pthread.h>
#include <math.h>

#define MAX_SIZE 4096

typedef double matrix[MAX_SIZE][MAX_SIZE];


int	N, n;		/* matrix size & threads*/
int	maxnum;		/* max number of element*/
char	*Init;		/* matrix init type	*/
int	PRINT;		/* print switch		*/
matrix	A;		/* matrix A		*/
double	b[MAX_SIZE];	/* vector b             */
double	y[MAX_SIZE];	/* vector y             */

// Structure to store the thread parameters
struct s1 {
int norm;
int id;
int block;
};

/* forward declarations */
void work();
void* calcRow(void *s);
void* calcBlocks(void *s);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);
pthread_barrier_t barr;  // declaration of barrier variable
int 
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;
 
    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    work();
    if (PRINT == 1)
	Print_Matrix();
}

void work() {
  int k;  //Normalization row
  printf("Computing Parallel via Pthreads. n=%d\n", n);

  
  for (k = 0; k < N; k++) {
       pthread_t threads[n];
	int i,j;
	struct s1* para = malloc(n * sizeof(struct s1));
	int block = (N-(k+1)+(n-1))/n; // calculation of ceil((N-(k+1))/n)
				     // but ceil() func gave compile time error 
				     // even when <math.h> is included. 

	pthread_barrier_init(&barr, NULL, n); // Initiating the barrier
	for (i = 0; i < n; i++) {
		para[i].norm = k;
		para[i].id = i;
		para[i].block=block;
		/* Creating n threads running calcRow function with argument &para[i] */
		pthread_create(&threads[i], NULL, calcRow, (void*) &para[i]);
	}
	
	/* Joining all threads */
	for (j = 0; j < n; j++)
		pthread_join(threads[j], NULL);
	
	y[k] = b[k] / A[k][k];
	A[k][k] = 1.0;
	
	//Distributing the work again but using the same parameters stucture for threads as last time
	for (i = 0; i < n; i++)
		pthread_create(&threads[i], NULL, calcBlocks, (void*) &para[i]);
	
	/* Joining all threads again */
	for (j = 0; j < n; j++)
		pthread_join(threads[j], NULL);
	

	
	/* Freeing the memory for array of struct s1 */
	free(para);
  }  
}





 /* A fucntion calculate every row when doing the division step   */
void* calcRow(void *s) {
    struct s1* myStruct = (struct s1*) s;
	int k = myStruct->norm;
	int id = myStruct->id;
	int block=myStruct->block;
	int j;
	int start=k+1 +id*block;
	int end=start+block;
	if(end > N) end=N;// If elements of row are not evenly distributable then last thread will get less work;
	
	for(j=start;j<end;j++) {
	    A[k][j] = A[k][j] / A[k][k];//division step
    	}
    	pthread_barrier_wait(&barr);// Sync all threads
	pthread_exit(0);
}

// Function to calculate blocks of entire rows during the elimination step
void* calcBlocks(void *s){
	struct s1* myStruct = (struct s1*) s;
	int k = myStruct->norm;
	int id = myStruct->id;
	int block=myStruct->block;
	int i,j;
	int start=k+1 +id*block;
	int end=start+block;
	if(end > N) end=N;// If the rows are not evenly distributable then last thread will get less rows;
	for (i = start; i < end; i++) {
	    for (j = k+1; j < N; j++)
		A[i][j] = A[i][j] - A[i][k]*A[k][j]; /* Elimination step */
	    b[i] = b[i] - A[i][k]*y[k];
	    A[i][k] = 0.0;
	}
	pthread_barrier_wait(&barr); // Sync all threads
	pthread_exit(0);
}

void
Init_Matrix()
{
    int i, j;
 
    printf("\nsize      = %dx%d ", N, N);
    printf("\nmaxnum    = %d \n", maxnum);
    printf("Init	  = %s \n", Init);
    printf("Initializing matrix...");
 
    if (strcmp(Init,"rand") == 0) {
	for (i = 0; i < N; i++){
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = (double)(rand() % maxnum) + 5.0;
		else
		    A[i][j] = (double)(rand() % maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(Init,"fast") == 0) {
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = 5.0;
		else
		    A[i][j] = 2.0;
	    }
	}
    }

    /* Initialize vectors b and y */
    for (i = 0; i < N; i++) {
	b[i] = 2.0;
	y[i] = 1.0;
    }

    printf("done \n\n");
    if (PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;
 
    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
	printf("[");
	for (j = 0; j < N; j++)
	    printf(" %5.2f,", A[i][j]);
	printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", y[j]);
    printf("]\n");
    printf("\n\n");
}

void 
Init_Default()
{
    N = 2048;
    Init = "rand";
    maxnum = 5.0;
    PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'N':
		--argc;
		n = atoi(*++argv);
		break;
		case 'n':
		--argc;
		N = atoi(*++argv);
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", N);
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		Init = *++argv;
		break;
	    case 'm':
		--argc;
		maxnum = atoi(*++argv);
		break;
	    case 'P':
		--argc;
		PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
