/*****************************************************
 *
 * S O R algorithm
 * ("Red-Black" solution to LaPlace approximation)
 *
 * Parallel version
 *
 *****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#define MAX_SIZE 4096
#define EVEN_TURN 0 /* shall we calculate the 'red' or the 'black' elements */
#define ODD_TURN  1
#define FROM_MASTER 1
#define FROM_WORKER_EVEN 20
#define FROM_WORKER 2
#define FROM_WORKER_ODD 21
MPI_Status status;

typedef double matrix[MAX_SIZE+2][MAX_SIZE+2]; /* (+2) - boundary elements */

 struct globmem {
    int		N;		/* matrix size		*/
    int		maxnum;		/* max number of element*/
    char	*Init;		/* matrix init type	*/
    double	difflimit;	/* stop condition	*/
    double	w;		/* relaxation factor	*/
    int		PRINT;		/* print switch		*/
    matrix	A;		/* matrix A		*/
} *glob;

/* forward declarations */
int work();
void Init_Matrix();
void Print_Matrix();
void Init_Default();
int Read_Options(int, char **);
int proc;
int 
main(int argc, char **argv)
{
    int timestart, timeend,iter; 
    
    
    glob = (struct globmem *) malloc(sizeof(struct globmem));
    Init_Default();		/* Init default values	*/
	Read_Options(argc,argv);	/* Read arguments	*/
	Init_Matrix();		/* Init the matrix	*/
	MPI_Init(&argc, &argv);
	timestart = MPI_Wtime();
        iter = work();
        timeend = MPI_Wtime();
        if (glob->PRINT == 1)
            Print_Matrix();
        printf("\nNumber of iterations = %d\n", iter);
        printf("Execution time on %2d nodes: %2d\n", proc, timeend-timestart);
        
        MPI_Finalize();
}
int 
work()
{
    int nproc, myrank, turn=EVEN_TURN,rows, offset, mtype, N, i, j, m, n, finished=0, iteration=0;
    double prevmax_even, prevmax_odd, maxi, sum, w, least; 
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	proc=nproc;
    if(myrank == 0)
    {//printf("master entered");
        prevmax_even = 0.0;
        prevmax_odd = 0.0;
	N = glob->N;
	w = glob->w;
	mtype = FROM_MASTER;
	rows = N/nproc;
	offset = rows;
	
		
		
	for(i=1;i<nproc;i++)
	{
		
		MPI_Send(&rows, 1, MPI_INT, i, mtype, MPI_COMM_WORLD);
		MPI_Send(&offset, 1, MPI_INT, i, mtype, MPI_COMM_WORLD);
		MPI_Send(&glob->A[offset][0], (rows+2)*(N+2), MPI_DOUBLE, i, mtype, MPI_COMM_WORLD);
      		offset += rows;
		//printf("\nfirst send mblock\n");
	}
		
	
	 while (!finished)
	{	 
            iteration++;
            if(turn==EVEN_TURN)
            {
	    	for (m = 1; m < rows+1; m++) 
		{
                    for (n = 1; n < N+1; n++) 
                    {
		    	if (((m + n) % 2) == 0)
                            glob->A[m][n] = (1 - w) * glob->A[m][n] + w * (glob->A[m-1][n] + glob->A[m+1][n] + glob->A[m][n-1] + glob->A[m][n+1]) / 4;
                    }
	    	}
	    	/* Calculate the maximum sum of the elements */
	    	maxi = -999999.0;
	    	for (m = 1; m < rows+1; m++) 
                {
                    sum = 0.0;
                    for (n = 1; n < N+1; n++)
		    	sum += glob->A[m][n];
                    if (sum > maxi)
		    	maxi = sum;
	    	}
	//	printf("%f", maxi);
		least=maxi;
		mtype = FROM_WORKER_EVEN;
		for(i=1;i<nproc;i++)
                {
                    //printf("before first mreceive");
                    MPI_Recv(&maxi, 1, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                    if(least<=maxi)
			least = maxi;
                //    printf("%f\n",maxi);
		}
		
                /* Compare the sum with the prev sum, i.e., check wether 
                * we are finished or not. */

	    	if (fabs(least - prevmax_even) <= glob->difflimit)
		    finished = 1;
		//printf("fin=%d\n",finished);

                if ((iteration%100) == 0)
		    printf("Iteration: %d, least = %f, prevmax_even = %f\n",iteration, least, prevmax_even);
		prevmax_even=least;
		turn = ODD_TURN;
            }
            else if(turn==ODD_TURN)
            {
		/* CALCULATE part B - odd elements*/
                for (m = 1; m < rows+1; m++) 
                {
                    for (n = 1; n < N+1; n++) 
                    {
                        if (((m + n) % 2) == 1)
                            glob->A[m][n] = (1 - w) * glob->A[m][n] + w * (glob->A[m-1][n] + glob->A[m+1][n] + glob->A[m][n-1] + glob->A[m][n+1]) / 4;
                    }
                }
                /* Calculate the maximum sum of the elements */
                maxi = -999999.0;
                for (m = 1; m < rows+1; m++) 
                {
                    sum = 0.0;
                    for (n = 1; n < N+1; n++)
                        sum += glob->A[m][n];
                    if (sum > maxi)
                        maxi = sum;
                }
		
                least = maxi;
                mtype = FROM_WORKER_ODD;
                for(i=1;i<nproc;i++)
                {
                    //printf("before second mreceive");
                    MPI_Recv(&maxi, 1, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                    if(least<=maxi)
                        least = maxi;
                    //printf("second mreceive\n");
                }
		
                /* Compare the sum with the prev sum, i.e., check wether 
                * we are finished or not. */

                if (fabs(least - prevmax_odd) <= glob->difflimit)
                    finished = 1;
	//	printf("fin=%d\n",finished);
                if ((iteration%100) == 0)
                    printf("Iteration: %d, least = %f, prevmax_odd = %f\n",iteration, least, prevmax_odd);
                prevmax_odd=least;
	        turn = EVEN_TURN;
            }
            
            if (iteration > 100000) 
            {
                /* exit if we don't converge fast enough */
                printf("Max number of iterations reached! Exit!\n");
                finished = 1;
            }
	 
            // data from workers
            mtype = FROM_WORKER;
            for(i=1;i<nproc;i++)
            {	
		MPI_Recv(&rows, 1, MPI_INT, i, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&offset, 1, MPI_INT, i, mtype, MPI_COMM_WORLD, &status);
		
		for(m=offset+1;m<offset+rows+1;m++)
                    MPI_Recv(&glob->A[m][0], N+2, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                    //printf("third mreceive");
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
	}//MPI_Finalize();
	MPI_Barrier(MPI_COMM_WORLD);
	
	return iteration;
		
    }	
    else // Worker part
    {	
	N = glob->N;
	w = glob->w;
	mtype = FROM_MASTER;
		
        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
      	MPI_Recv(&offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&glob->A[offset][0], (rows+2)*(N+2), MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
        //printf("first wreceive");
				
	//calculation by worker
        while (!finished)
	{
            iteration++;
            /* CALCULATE part A - even elements */
            if (turn == EVEN_TURN)
            {
                for (m = (offset+1); m < (offset+rows+2); m++) 
		{
                    for (n = 1; n < N+1; n++) 
                    {
			if (((m + n) % 2) == 0)
                            glob->A[m][n] = (1 - w) * glob->A[m][n]	+ w * (glob->A[m-1][n] + glob->A[m+1][n] + glob->A[m][n-1] + glob->A[m][n+1]) / 4;		
                    }
		}
                /* Calculate the maximum sum of the elements */
                maxi = -999999.0;
                for (m = (offset+1); m < (offset+rows+2); m++)
                {
                    sum = 0.0;
                    for (n = 1; n < N+1; n++)
                        sum += glob->A[m][n];
                    if (sum > maxi)
                        maxi = sum;
                } 
	    
                /* Compare the sum with the prev sum, i.e., check wether 
                * we are finished or not. */
                mtype = FROM_WORKER_EVEN;
		//printf("before first wsend");
		MPI_Send(&maxi, 1, MPI_DOUBLE, 0, mtype,MPI_COMM_WORLD);
		//printf("first wsend");
		
		
		
	
                    turn=ODD_TURN;
            } 
            else if (turn == ODD_TURN) 
            {
		for (m = (offset+1); m < (offset+rows+2); m++) 
		{
                    for (n = 1; n < N+1; n++) 
                    {
			if (((m + n) % 2) == 1)
                            glob->A[m][n] = (1 - w) * glob->A[m][n]	+ w * (glob->A[m-1][n] + glob->A[m+1][n] + glob->A[m][n-1] + glob->A[m][n+1]) / 4;		
                    }
		}
                /* Calculate the maximum sum of the elements */
                maxi = -999999.0;
                for (m = (offset+1); m < (offset+rows+2); m++)
                {
                    sum = 0.0;
                    for (n = 1; n < N+1; n++)
                        sum += glob->A[m][n];
                    if (sum > maxi)
                        maxi = sum;
		} 
	    
                /* Compare the sum with the prev sum, i.e., check wether 
                * we are finished or not. */
                mtype = FROM_WORKER_ODD;
		//printf("before second wsend");
		MPI_Send(&maxi, 1, MPI_DOUBLE, 0, mtype,MPI_COMM_WORLD);
		//printf("second wsend");
		
		
		
		
			turn=EVEN_TURN;
            }
            mtype = FROM_WORKER;
		MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
            
            for(m=offset+1;m<offset+rows+1;m++)
		MPI_Send(&glob->A[m][0], N+2, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
		//printf("third wsend");
            
            MPI_Barrier(MPI_COMM_WORLD);
	}
        MPI_Barrier(MPI_COMM_WORLD);
    }	
	//MPI_Finalize();
}

void
Init_Matrix()
{
    int i, j, N, dmmy;
 
    N = glob->N;
    printf("\nsize      = %dx%d ",N,N);
    printf("\nmaxnum    = %d \n",glob->maxnum);
    printf("difflimit = %.7lf \n",glob->difflimit);
    printf("Init	  = %s \n",glob->Init);
    printf("w	  = %f \n\n",glob->w);
    printf("Initializing matrix...");
 
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i < glob->N+2; i++) {
	for (j = 0; j < glob->N+2; j++) {
	    glob->A[i][j] = 0.0;
	}
    }
    if (strcmp(glob->Init,"count") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		glob->A[i][j] = (double)i/2;
	    }
	}
    }
    if (strcmp(glob->Init,"rand") == 0) {
	for (i = 1; i < N+1; i++){
	    for (j = 1; j < N+1; j++) {
		glob->A[i][j] = (rand() % glob->maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(glob->Init,"fast") == 0) {
	for (i = 1; i < N+1; i++){
	    dmmy++;
	    for (j = 1; j < N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    glob->A[i][j] = 1.0;
		else
		    glob->A[i][j] = 5.0;
	    }
	}
    }

    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    glob->A[0][0] = glob->A[1][1];
    glob->A[0][N+1] = glob->A[1][N];
    glob->A[N+1][0] = glob->A[N][1];
    glob->A[N+1][N+1] = glob->A[N][N];
    /* fix the top and bottom rows */
    for (i = 1; i < N+1; i++) {
	glob->A[0][i] = glob->A[1][i];
	glob->A[N+1][i] = glob->A[N][i];
    }
    /* fix the left and right columns */
    for (i = 1; i < N+1; i++) {
	glob->A[i][0] = glob->A[i][1];
	glob->A[i][N+1] = glob->A[i][N];
    }

    
    if (glob->PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix()
{
    int i, j, N;
 
    N = glob->N;
    for (i=0; i<N+2 ;i++){
	for (j=0; j<N+2 ;j++){
	    printf(" %f",glob->A[i][j]);
	}
	printf("\n");
    }
    printf("\n\n");
}

void 
Init_Default()
{
    glob->N = 2048;
    glob->difflimit = 0.00001*glob->N;
    glob->Init = "rand";
    glob->maxnum = 15.0;
    glob->w = 0.5;
    glob->PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		glob->N = atoi(*++argv);
		glob->difflimit = 0.00001*glob->N;
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-d difflimit] 0.1-0.000001 \n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand/count \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		printf("           [-w relaxation_factor] 1.0-0.1 \n\n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", glob->N);
		printf("\n          difflimit = 0.0001 ");
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          w         = 0.5 \n");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		glob->Init = *++argv;
		break;
	    case 'm':
		--argc;
		glob->maxnum = atoi(*++argv);
		break;
	    case 'd':
		--argc;
		glob->difflimit = atof(*++argv);
		break;
	    case 'w':
		--argc;
		glob->w = atof(*++argv);
		break;
	    case 'P':
		--argc;
		glob->PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}

