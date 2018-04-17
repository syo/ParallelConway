/***************************************************************************/
/* Asssignment 4/5              ********************************************/
/* Connor Hadley                **(*****************************************/
/* Zach Biletch                 **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>
#include<pthread.h>

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void *threadcall(void *void_ptr) {
    printf("called thread\n");
}

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int gridsize = 16384;
    int mpi_myrank;
    int mpi_commsize;
    double start, end;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
    // Start timing
    if (mpi_myrank == 0) {
        start = MPI_Wtime();
    }

    // Allocate rank's rows and ghost rows for the boundary
    char** rows = calloc(gridsize / mpi_commsize, sizeof(char));
    char* rowbefore = calloc(gridsize, sizeof(char));
    char* rowafter = calloc(gridsize, sizeof(char));
    for (int i = 0; i < gridsize; i++) {
        rows[i] = calloc(gridsize, sizeof(char));
    }
    
//XXX
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code

    /* THREAD STUFF */
    /*
        Will need to create an array of pthread_t for each MPI rank
        Should be relatively straightforward
        Can use this basic framework to operate, need to specialize the function to modify global variable for shared memory maybe?
    */
    int val;
    pthread_t somethread;
    // create a thread to run function threadcall
    if (pthread_create(&somethread, NULL, threadcall, &val)) {
        fprintf(stderr, "Error creating thread\n");
        return 1;
    }

    //do whatever main thread needs to do
    printf("from thread 1\n");

    //wait for the thread to finish and join it back in
    if(pthread_join(somethread, NULL)) {
        fprintf(stderr, "Error joining thread\n");
        return 2;
    }

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
