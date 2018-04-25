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

/*#define BGQ 1
#ifdef BGQ
#include <hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif*/
#define ALIVE '1'
#define DEAD  '0'
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
pthread_barrier_t mpi_io;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void *threadcall(void *val_ptr) {
    //printf("called thread\n");
    int base_thread = pthread_mutex_trylock(&threadzero); //use a mutex to make sure only one thread does mpi things
    //if (base_thread == 0)
    //    printf("this should show up 4 times\n"); 
    // Overall for loop
    for (int i = 0; i < numticks; i++) {
        //printf("TICK %d\n", i);
        // MPI send and receive
        if (base_thread == 0) {
            //printf("started send and received at TICK %d\n", i);
            // Recieve the row before and after
            MPI_Request recva, recvb, senda, sendb;
            int before = (mpi_myrank + mpi_commsize - 1) % mpi_commsize;
            int after = (mpi_myrank + 1) % mpi_commsize;
            MPI_Irecv(rowbefore, gridsize, MPI_CHAR, before, 123, MPI_COMM_WORLD, &recvb);
            MPI_Irecv(rowafter, gridsize, MPI_CHAR, after, 321, MPI_COMM_WORLD, &recva);

            // Send the first and last rows to neighbors
            MPI_Isend(rows[0], gridsize, MPI_CHAR, before, 321, MPI_COMM_WORLD, &sendb);
            MPI_Isend(rows[rowsperrank-1], gridsize, MPI_CHAR, after, 123, MPI_COMM_WORLD, &senda);

            // Wait for send and receive to finish
            MPI_Status status;
            MPI_Wait(&recvb, &status);
            MPI_Wait(&recva, &status);
            MPI_Wait(&sendb, &status);
            MPI_Wait(&senda, &status);
        }

        pthread_barrier_wait(&mpi_io);
        //printf("did mpi io\n");
        // process each row belonging to this pthread
        // should pass in what belongs to it as an argument to threadcall
        // use modulus to calculate positions
        int cur_row = *((int *) val_ptr);
        int j;
        for (j=0; j<rowsperthread; j++) {//iterate through each row in the thread
            int lock_status = pthread_mutex_trylock(&gridlock);//if we can get the mutex, process a row
            while (lock_status == EBUSY) {
                // if mutex is locked just wait for it
                lock_status = pthread_mutex_trylock(&gridlock);
            }
            //printf("LOCKED MUTEX\n");
            // write new values for the rows
            // need to possibly calculate beforehand so updates don't change outcomes of other threads?

            //calculate each new value in the row
            int k;
            for (k=0; k < gridsize; k++) {//for item in row
                //printf("infinite? %d\n", k);
                int life_status = rows[cur_row + j][k] - '0';
                int neighbors = 0; // XXX need to calc this
                if (k > 0) {//check to the left
                    neighbors += rows[cur_row + j][k - 1] - '0';
                }
                if (k < gridsize - 1) {//check to the right
                    neighbors += rows[cur_row + j][k + 1] - '0';
                }
                if (j > 0) { //if this is not the first row of the thread
                    //check above
                    neighbors += rows[cur_row + j - 1][k - 1] - '0';
                    neighbors += rows[cur_row + j - 1][k] - '0';
                    neighbors += rows[cur_row + j - 1][k + 1] - '0';
                } else { //if this is the first row of the thread
                    if (cur_row == 0) {//check if its the first thread of mpi rank
                        //ask for the ghost row info
                        neighbors += rowbefore[k - 1] - '0'; 
                        neighbors += rowbefore[k] - '0'; 
                        neighbors += rowbefore[k + 1] - '0'; 
                    } else {
                        //ask for the thread above's row XXX may be doing this wrong
                        neighbors += rows[cur_row + j - 1][k - 1] - '0';
                        neighbors += rows[cur_row + j - 1][k] - '0';
                        neighbors += rows[cur_row + j - 1][k + 1] - '0';
                    }
                }
                if (j < rowsperthread - 1) { //if this is not the last row of the thread
                    //check below 
                    neighbors += rows[cur_row + j + 1][k - 1] - '0';
                    neighbors += rows[cur_row + j + 1][k] - '0';
                    neighbors += rows[cur_row + j + 1][k + 1] - '0';
                } else { //if this is the last row of the thread
                    //printf("last row of thread\n");
                    if (cur_row / rowsperthread + 1 == THREADS) {//check if its the last thread of mpi rank
                        //ask for the ghost row info
                        //printf("last thread of rank\n");
                        neighbors += rowafter[k - 1] - '0'; 
                        neighbors += rowafter[k] - '0'; 
                        neighbors += rowafter[k + 1] - '0'; 
                    } else {
                        //ask for the thread below's row XXX may be doing this wrong
                        neighbors += rows[cur_row + j + 1][k - 1] - '0';
                        neighbors += rows[cur_row + j + 1][k] - '0';
                        neighbors += rows[cur_row + j + 1][k + 1] - '0';
                    }
                }
                if (life_status == 1) {
                    //if alive
                    if (neighbors < 2 || neighbors > 3) {// if not 2 or 3 neighbors
                        rows[cur_row + j][k] = '0'; //kill it
                    }
                } else {
                    // if dead
                    if (neighbors == 3) { //if exactly 3 neighbors
                        rows[cur_row + j][k] = '1'; //bring it to life
                    }
                }
                //if (i == 0) 
                //    printf("updated value from -%c- to -%c-\n", life_status, rows[cur_row + j][k]);
            }
            //printf("UNLOCKED MUTEX\n");
            pthread_mutex_unlock(&gridlock);
        }
        //printf("updated rows in thread\n");
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
    //printf("running program\n");
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    rowsperrank = gridsize / mpi_commsize;
    rowsperthread = rowsperrank / THREADS;

    //init grid mutex
    pthread_mutex_init(&gridlock, NULL);

    //init threadzero mutex
    pthread_mutex_init(&threadzero, NULL);

    //init pthread barrier
    pthread_barrier_init(&mpi_io, NULL, THREADS);

    // Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
    // Start timing
    if (mpi_myrank == 0) {
        printf("starting timer\n");
        start = MPI_Wtime();
    }
    if (mpi_myrank == 0) {
        printf("timer started in rank 0\n");
    }

    // Allocate rank's rows and ghost rows for the boundary
    rows = calloc(rowsperrank, sizeof(char*));
    rowbefore = calloc(gridsize, sizeof(char));
    rowafter = calloc(gridsize, sizeof(char));
    for (int i = 0; i < gridsize; i++) {
        rows[i] = calloc(gridsize, sizeof(char));
        //rows[i][gridsize] = '\0';
    }

    //rows[gridsize] = '\0';
    //rowafter[gridsize] = '\0';
    //rowbefore[gridsize] = '\0';

    if (mpi_myrank == 0) {
        printf("allocated in rank 0\n");
    }

    //printf("rows allocated, initializing\n");
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
    if (mpi_myrank == 0) {
        printf("initialized in rank 0\n");
    }
    //printf("initialized\n");
//XXX
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    //printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	//   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier(MPI_COMM_WORLD);
    

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

    if (mpi_myrank == 0) {
        printf("creating threads in rank 0\n");
    }
    // create the threads calling the threadcall function
    int i;
    for (i=0; i<num_threads; i++) {
        int val = currow; //pass it the current row as the beginning of this thread's rows
        pthread_create(&p_threads[i], &attr, threadcall, (void *) &val);
        currow += rowsperthread; //increment current row by number of rows per thread
    }

    //do whatever main thread needs to do
    //printf("this is the main thread\n");

    if (mpi_myrank == 0) {
        printf("joining threads in rank 0\n");
    }
    //wait for threads to finish and join them back in
    for (i=0; i<num_threads; i++) {
        pthread_join(p_threads[i], NULL);
    }

    // Sync after threads end
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_myrank == 0) {
        // End of compute timing
        end = MPI_Wtime();
        comptime = end - start;

        // I/O timing
        start = MPI_Wtime();
    }
    
    if (mpi_myrank == 0) {
        printf("performing io in rank 0\n");
    }
    // Perform I/O
    MPI_File outfile;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile);
    for (int i = 0; i < rowsperrank; i++) {
        //rowtoout(i);
        MPI_File_write_at(outfile, (mpi_myrank * rowsperrank + i) * (gridsize + 1), rows[i], gridsize, MPI_CHAR, &status);
        // Add newline at the end of each line
        char newline = '\n';
        MPI_File_write_at(outfile, (mpi_myrank * rowsperrank + i + 1) * (gridsize + 1) - 1, &newline, 1, MPI_CHAR, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&outfile);

    if (mpi_myrank == 0) {
        // End of I/O timing
        end = MPI_Wtime();
        iotime = end - start;

        printf("Compute time: %lfs, I/O time: %lfs\n", comptime, iotime);
    }

    // Output heatmap data to a file
    char heatmap = 1;
    if (heatmap) {
        MPI_File heatfilel;
        MPI_File heatfiler;
        //if (mpi_myrank < mpi_commsize / 2) {
            printf("opening files 1 and 2\n");
            MPI_File_open(MPI_COMM_WORLD, "heatmap1.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &heatfilel);
            MPI_File_open(MPI_COMM_WORLD, "heatmap2.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &heatfiler);
            printf("opened files %p and %p\n", heatfilel, heatfiler);
            /*
        } else {
            printf("opening files 3 and 4\n");
            MPI_File_open(MPI_COMM_WORLD, "heatmap3.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &heatfilel);
            MPI_File_open(MPI_COMM_WORLD, "heatmap4.dat", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &heatfiler);
            printf("opened files %p and %p\n", heatfilel, heatfiler);
        }
        */


        int heatgridsize = 1024;
        int heatrowsperrank = 1024 / mpi_commsize;
        int linesperrank = heatgridsize * heatrowsperrank / 2;
        char* line = calloc(32, sizeof(char));
        for (int i = 0; i < heatrowsperrank; i++) {
            for (int j = 0; j < heatgridsize * 2; j++) {
                int sum = 0;
                for (int x = 0; x < 8; x++) {
                    for (int y = 0; y < 16; y++) {
                        sum += rows[i*8+x][j*16+y] - '0';
                    }
                }
                int numinline = sprintf(line, "%4d %4d %3d\n", (heatrowsperrank * mpi_myrank) + i, j, sum);
                if (j < heatgridsize) {
                    MPI_File_write_at(heatfilel, ((linesperrank * mpi_myrank) + (i * heatgridsize) + j) * numinline, line, numinline, MPI_CHAR, &status);
                } else {
                    MPI_File_write_at(heatfiler, ((linesperrank * mpi_myrank) + (i * heatgridsize) + (j - heatgridsize)) * numinline, line, numinline, MPI_CHAR, &status);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_close(&heatfilel);
        MPI_File_close(&heatfiler);
        free(line);
    }


    printf("rank %d rows at %p points to: %p\n", mpi_myrank, (void*)&rows, (void*)rows);
    // Clean up environment
    for (int i = 0; i < rowsperrank; i++) {
        //printf("trying to free row %d\n", i);
        free(rows[i]);
    }
    printf("freed each row\n");
    free(rows);
    printf("freed rows\n");
    free(rowbefore);
    printf("freed rowbefore\n");
    free(rowafter);
    printf("freed rowafter\n");
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
