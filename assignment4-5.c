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
#define THREADS 4

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int gridsize = 16384;
int numticks = 128;
int mpi_myrank;
int mpi_commsize;
int rowsperrank;
int rowsperthread;
char** rows;
char* rowbefore;
char* rowafter;
pthread_mutex_t gridlock;
pthread_mutex_t threadzero;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void *threadcall(void *val_ptr) {
    printf("called thread\n");
    int base_thread = pthread_mutex_trylock(&threadzero); //use a mutex to make sure only one thread does mpi things

    // Overall for loop
    for (int i = 0; i < numticks; i++) {
        // MPI send and receive
        if (base_thread == 0) {
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

        // process each row belonging to this pthread
        // should pass in what belongs to it as an argument to threadcall
        // use modulus to calculate positions
        int cur_row = *((int *) val_ptr);
        int j;
        for (j=0; j<rowsperthread; j++) {//iterate through each row in the thread
            int lock_status = pthread_mutex_trylock(&gridlock);
            if (lock_status == EBUSY) {
                // if mutex is locked just try this row again until it isn't
                j -= 1;
                continue;
            }
            else {
                // write new values for the rows
                // need to possibly calculate beforehand so updates don't change outcomes of other threads?

                //calculate each new value in the row
                int k;
                for (k=0; k < gridsize; k++) {//for item in row
                    int life_status = rows[cur_row + j][k];
                    int neighbors = 0; // XXX need to calc this
                    if (k > 0) {//check to the left
                        neighbors += rows[cur_row + j][k - 1];
                    }
                    if (k < gridsize - 1) {//check to the right
                        neighbors += rows[cur_row + j][k + 1];
                    }
                    if (j > 0) { //if this is not the first row of the thread
                        //check above
                        neighbors += rows[cur_row + j - 1][k - 1];
                        neighbors += rows[cur_row + j - 1][k];
                        neighbors += rows[cur_row + j - 1][k + 1];
                    } else { //if this is the first row of the thread
                        if (cur_row == 0) {//check if its the first thread of mpi rank
                            //ask for the ghost row info
                            neighbors += rowbefore[k - 1]; 
                            neighbors += rowbefore[k]; 
                            neighbors += rowbefore[k + 1]; 
                        } else {
                            //ask for the thread above's row XXX may be doing this wrong
                            neighbors += rows[cur_row + j - 1][k - 1];
                            neighbors += rows[cur_row + j - 1][k];
                            neighbors += rows[cur_row + j - 1][k + 1];
                        }
                    }
                    if (j < rowsperthread) { //if this is not the last row of the thread
                        //check below 
                        neighbors += rows[cur_row + j + 1][k - 1];
                        neighbors += rows[cur_row + j + 1][k];
                        neighbors += rows[cur_row + j + 1][k + 1];
                    } else { //if this is the last row of the thread
                        if (cur_row / THREADS + 1 == rowsperthread) {//check if its the last thread of mpi rank
                            //ask for the ghost row info
                            neighbors += rowafter[k - 1]; 
                            neighbors += rowafter[k]; 
                            neighbors += rowafter[k + 1]; 
                        } else {
                            //ask for the thread below's row XXX may be doing this wrong
                            neighbors += rows[cur_row + j + 1][k - 1];
                            neighbors += rows[cur_row + j + 1][k];
                            neighbors += rows[cur_row + j + 1][k + 1];
                        }
                    }
                    if (life_status) {
                        //if alive
                        if (neighbors < 2 || neighbors > 3) {// if not 2 or 3 neighbors
                            rows[cur_row + j][k] = 0; //kill it
                        }
                    } else {
                        // if dead
                        if (neighbors == 3) { //if exactly 3 neighbors
                            rows[cur_row + j][k] = 1; //bring it to life
                        }
                    }
                }
                pthread_mutex_unlock(&gridlock);
            }
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
    printf("running program\n");
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    rowsperrank = gridsize / mpi_commsize;
    rowsperthread = rowsperrank / THREADS;

    //init grid mutex
    pthread_mutex_init(&gridlock, NULL);

    //init threadzero mutex
    pthread_mutex_init(&threadzero, NULL);

    // Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
    // Start timing
    if (mpi_myrank == 0) {
        printf("starting timer\n");
        start = MPI_Wtime();
    }

    // Allocate rank's rows and ghost rows for the boundary
    rows = calloc(rowsperrank, sizeof(char*));
    rowbefore = calloc(gridsize, sizeof(char));
    rowafter = calloc(gridsize, sizeof(char));
    for (int i = 0; i < gridsize; i++) {
        rows[i] = calloc(gridsize, sizeof(char));
    }

    printf("rows allocated, initializing\n");
    // Randomly initialize universe
    for (int i = 0; i < rowsperrank; i++) {
        for (int j = 0; j < gridsize; j++) {
            if (GenVal(i) < 0.5) {
                rows[i][j] = DEAD;
            } else {
                rows[i][j] = ALIVE;
            }
        }
    }
    
    printf("initialized\n");
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

    // each rank intiializes an array of pthreads
    pthread_t p_threads[THREADS];
    int currow = 0; //start with the first row of the rank
    int num_threads = sizeof(p_threads) / sizeof(pthread_t);

    //init attr
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    // create the threads calling the threadcall function
    int i;
    for (i=0; i<num_threads; i++) {
        int val = currow; //pass it the current row as the beginning of this thread's rows
        pthread_create(&p_threads[i], &attr, threadcall, (void *) &val);
        currow += rowsperthread; //increment current row by number of rows per thread
    }

    //do whatever main thread needs to do
    //printf("this is the main thread\n");

    //wait for threads to finish and join them back in
    for (i=0; i<num_threads; i++) {
        pthread_join(p_threads[i], NULL);
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
        MPI_File_write_at(outfile, (mpi_myrank * rowsperrank + i + 1) * (gridsize + 1) - 1, &newline, 1, MPI_CHAR, &status);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_File_close(&outfile);

    if (mpi_myrank == 0) {
        // End of I/O timing
        end = MPI_Wtime();
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
