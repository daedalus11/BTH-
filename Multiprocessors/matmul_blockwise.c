
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 8	/* assumption: SIZE a multiple of number of nodes */
#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */
#define DEBUG 1		/* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];
static double bt[SIZE][SIZE];

static void
init_matrix(void)
{
    int i, j;
// intializing left and top half of A and B respectively as 1.0 and the other halfs as 2.0
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
		if(i < (SIZE/2))
	    		a[i][j] = 1.0;
		else 
			a[i][j] = 2.0;
		if(j < (SIZE/2))
	    		b[i][j] = 1.0;
		else
			b[i][j] = 2.0;
        }
for (i = 0; i < SIZE; i++) 		// transposing matrix B for easier distribution of columns as rows
        for (j = 0; j < SIZE; j++)
		bt[i][j]=b[j][i];
if(DEBUG){
    printf("*******MATRIX A**************\n");
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
	    printf(" %7.2f", a[i][j]);
	printf("\n");
    }
printf("*******MATRIX B**************\n");
for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
	    printf(" %7.2f", b[i][j]);
	printf("\n");
    }
	}

}

static void
print_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
	    printf(" %7.2f", c[i][j]);
	printf("\n");
    }
}

int
main(int argc, char **argv)
{
    int myrank, nproc;
    int rows, cols; /* amount of work per node (rows per worker) */
    int mtype; /* message type: send/recv between master and workers */
    int dest, src, roffset, coffset;
    double start_time, end_time;
    int i, j, k, blockrows, blockcolumns;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) {
	/* Master task */

	/* Initialization */
	printf("\nSIZE = %d, number of nodes = %d\n", SIZE, nproc);
	init_matrix();
	start_time = MPI_Wtime();

	/* Send part of matrix a and matrix b to workers */
	if (nproc==1)
        {blockrows=1; blockcolumns=1;} else if (nproc==2){blockrows=2; blockcolumns=1;}
        else if (nproc==4)
        {blockrows=2; blockcolumns=2;}
        else if (nproc==8)
        {blockrows=4; blockcolumns=2;}
        else { printf("processors invalid"); exit(0);}
		rows=SIZE/blockrows;
        cols= SIZE/blockcolumns;

	mtype = FROM_MASTER;
	roffset = 0;
	coffset = cols;
	
	for (dest = 1; dest < nproc; dest++) {
	    //if (DEBUG)
		//printf("   \nsending %d rows to task %d\n",rows,dest);
		if(nproc==2) roffset+=rows;
	    MPI_Send(&roffset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);	   
		MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
		MPI_Send(&cols, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
	    MPI_Send(&a[roffset][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
	   
		if(dest%2==0)
		{	 coffset = 0;
			 MPI_Send(&coffset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                         MPI_Send(&bt[coffset][0], cols*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);

		}
		else
		{	coffset = cols;if(nproc==2) coffset=0;
			 MPI_Send(&coffset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			 MPI_Send(&bt[coffset][0], cols*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
			roffset += rows;
		}

	   
	}
	/* let master do its part of the work */
	for (i = 0; i < rows; i++) {
	    for (j = 0; j < cols; j++) {
		c[i][j] = 0.0;
		for (k = 0; k < SIZE; k++)
		    c[i][j] = c[i][j] + a[i][k] * b[k][j];
	    }
	}

	/* collect the results from all the workers */
	mtype = FROM_WORKER;
	for (src = 1; src < nproc; src++) {
	    MPI_Recv(&roffset, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&cols, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&coffset, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
		for(i=roffset;i<roffset+rows;i++)
	    MPI_Recv(&c[i][coffset], cols, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
		//if (DEBUG)
		//printf("   \nrecvd %d rows from task %d, offset = %d\n",
		  //     rows, src, roffset);
	}

	end_time = MPI_Wtime();
	if (DEBUG){
	    /* Prints the resulting matrix c */
	printf("---------------MATRIX C ______________________\n");
	    print_matrix();}
	printf("\nExecution time on %2d nodes: %f\n", nproc, end_time-start_time);

    } else {
	/* Worker tasks */

	/* Receive data from master */
	mtype = FROM_MASTER;
	MPI_Recv(&roffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&a[roffset][0], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&coffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&bt[coffset][0], SIZE*cols, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
	//if (DEBUG)
	//    printf ("\nRank=%d, offset=%d, row =%d, a[offset][0]=%e, b[0][0]=%e\n",
	//	    myrank, roffset, rows, a[roffset][0], bt[0][0]);

	/* do the workers part of the calculation */
	
	for(i=0; i<SIZE; i++){ printf("\n");
                for(j=coffset; j<coffset+cols; j++)
			{b[i][j]=bt[j][i];}}
	for (i=roffset; i<roffset+rows; i++)
	    for (j=coffset; j<coffset+cols; j++) {
		c[i][j] = 0.0;
		for (k=0; k<SIZE; k++)
		    c[i][j] = c[i][j] + a[i][k] * b[k][j];
	    }
	//if (DEBUG)
	 //   printf ("\nRank=%d, offset=%d, row =%d, c[offset][0]=%e\n",
	//	    myrank, roffset, rows, a[roffset][0]);

	/* send the results to the master */
	mtype = FROM_WORKER;
	MPI_Send(&roffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	MPI_Send(&cols, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	MPI_Send(&coffset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	for(i=roffset;i<roffset+rows;i++)
	MPI_Send(&c[i][coffset], cols, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);	
   }

    MPI_Finalize();
    return 0;
}


