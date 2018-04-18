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
#define THRESHOLD 25

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int gridsize = 16384;
int numticks = 128;
int mpi_myrank;
int mpi_commsize;
int rowsperrank;
char** rows;
char* rowbefore;
char* rowafter;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void *threadcall(void *void_ptr) {
    printf("called thread\n");
    // Overall for loop
    for (int i = 0; i < numticks; i++) {
        // MPI send and receive
        //XXX not sure how to get "thread 0" -- pthread_self() is a thing though
        if (thread is 0) {
            // Recieve the row before and after
            MPI_Request recv;
            MPI_Request send;
            int before = (mpi_myrank + gridsize - 1) % gridsize;
            int after = (mpi_myrank + 1) % gridsize;
            MPI_Irecv(rowbefore, gridsize, MPI_CHAR, before, 123, MPI_COMM_WORLD, &recv);
            MPI_Irecv(rowafter, gridsize, MPI_CHAR, after, 321, MPI_COMM_WORLD, &recv);

            // Send the first and last rows to neighbors
            MPI_Isend(rows[0], gridsize, MPI_CHAR, before, 321, MPI_COMM_WORLD, &send);
            MPI_Irecv(rows[rowsperrank-1], gridsize, MPI_CHAR, after, 123, MPI_COMM_WORLD, &send);

            // Wait for send and receive to finish
            MPI_Status status;
            MPI_Wait(&recv, &status);
            MPI_Wait(&send, &status);
        }
    }
    pthread_exit(0);
}

void rowtoout(int n) {
    for (int i = 0; i < gridsize; i++) {
        rows[n][i] += '0';
    }
}

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    double start, end, comptime, iotime;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    rowsperrank = gridsize / mpi_commsize;
    
    // Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
    // Start timing
    if (mpi_myrank == 0) {
        start = MPI_Wtime();
    }

    // Allocate rank's rows and ghost rows for the boundary
    rows = calloc(rewsperrank, sizeof(char));
    rowbefore = calloc(gridsize, sizeof(char));
    rowafter = calloc(gridsize, sizeof(char));
    for (int i = 0; i < gridsize; i++) {
        rows[i] = calloc(gridsize, sizeof(char));
    }

    // Randomly initialize universe
    for (int i = 0; i < rowsperrank; i++) {
        for (int j = 0; j < gridsize; j++) {
            // XXX not sure how to initialize the grid
            if (GenVal(i) < 0.5) {
                rows[i][j] = DEAD;
            } else {
                rows[i][j] = ALIVE;
            }
        }
    }
    
//XXX
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    

    /* THREAD STUFF */
    /*
        Will need to create an array of pthread_t for each MPI rank
        Should be relatively straightforward
        Can use this basic framework to operate, need to specialize the function to modify global variable for shared memory maybe?
    */
    pthread_t p_threads[4];
    int val;
    int num_threads = sizeof(p_threads) / sizeof(pthread_t);

    //init attr
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    // create however many threads calling necessary function
    int i;
    for (i=0; i<num_threads; i++) {
        pthread_create(&p_threads[i], &attr, threadcall, &val);
    }

    //do whatever main thread needs to do
    printf("from thread 1\n");

    //wait for threads to finish and join them back in
    for (i=0; i<num_threads; i++) {
        pthread_join(&p_threads[i], &attr, threadcall, &val);
    }

    // Sync after threads end
    MPI_Barrier( MPI_COMM_WORLD );

    if (mpi_myrank == 0) {
        // End of compute timing
        end = MPI_Wtime();
        comptime = end - start;

        // I/O timing
        start = MPI_Wtime();
    }
    
    // Perform I/O
    MPI_File outfile;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile);
    for (int i = 0; i < rowsperrank; i++) {
        rowtoout(i);
        MPI_File_write_at(outfile, (mpi_myrank * rowsperrank + i) * (gridsize + 1), rows[i], gridsize, MPI_CHAR, &status);
        // Add newline at the end of each line
        char newline = '\n';
        MPI_File_write_at(outfile, (mpi_myrank * rowsperrank + i + 1) * (gridsize + 1) - 1, &newline, MPI_CHAR, &status);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_File_close(outfile);

    if (mpi_myrank == 0) {
        // End of I/O timing
        end = MPI_Wtime;
        iotime = end - start;

        printf("Compute time: %lfs, I/O time: %lfs\n", comptime, iotime);
    }


    // Clean up environment
    free(rows);
    free(rowbefore);
    free(rowafter);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
