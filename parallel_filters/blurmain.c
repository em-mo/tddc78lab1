#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

void calcDispls(int xsize, int ysize, int numProc, int *displacements, int *sendCounts);
void constructPixelDataType(MPI_Datatype* pixelDataType);

int main (int argc, char ** argv) 
{
    int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    pixel target[MAX_PIXELS];
    struct timespec stime, etime;
    #define MAX_RAD 1000


    /* MPI INIT*/
    int ierr = MPI_Init(&argc, &argv);
    int myId, numberProc;
    
    int* displacements;
    int* sendCounts;
    MPI_Datatype pixelType;

    double stime, etime;

    MPI_Comm_size(MPI_COMM_WORLD, &numberProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    double w[MAX_RAD];

    /* Take care of the arguments */

    if (argc != 4) 
    {
        fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) 
    {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        MPI_Finalize();
        exit(1);
    }

    /* read file */
    if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0){
        MPI_Finalize();
        exit(1);
    }

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        MPI_Finalize();
        exit(1);
    }

    printf("Has read the image, generating coefficients\n");;

    if(myId = 0)
        stime = MPI_Wtime();
    
    /* filter */
    get_gauss_weights(radius, w);

    printf("Calling filter\n");

    constructPixelDataType(pixelType);

    displacements = (int*) malloc(numberProc * sizeof(int));
    sendCounts = (int*) malloc(numberProc * sizeof(int));

    if(myId == 0)
        target = (pixel*) malloc(MAX_PIXELS * sizeof(pixel));

    ystart = displacements[myId] / xsize;
    ystop = ystart + sendCounts[myId] / xsize;
    blurfilter(xsize, ysize, src, radius, w, ystart, ystop);

    MPI_Gatherv(&src[displacemts[myId]], sendCounts[myId], pixelType, target, sendCounts, displacements, pixelType, 0, MPI_COMM_WORLD);

    if(myId == 0)
    {
        etime = MPI_Wtime();
        printf("Total time: %.6f", etime - stime);

        printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
            1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

        /* write result */
        printf("Writing output file\n");

        if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0){
            MPI_Finalize();
            exit(1);
        }
    }
    MPI_Finalize();
    return(0);
}

void calcDispls(int xsize, int ysize, int numProc, int myId, int *displacements, int *sendCounts)
{
    int currentDisplacement = 0;
    int sendCount, i, rows, restRows;
    rows = ysize / numProc;
    restRows = ysize % numProc;

    for (i = 0; i < numProc; i++) 
    {
        displacements[i] = currentDisplacement;

        sendCount = restRows > 0 ? rows + 1 : rows;
        sendCount *= xsize;
        sendCounts[i] = sendCount;

        currentDisplacement += sendCount;

        restRows--;
    }
}

void constructPixelDataType(MPI_Datatype* pixelDataType) 
{
    MPI_Type_contiguous(3, MPI_CHAR, &pixelType);
    MPI_Type_commit(&pixelType);
}
