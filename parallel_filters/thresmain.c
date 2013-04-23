#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "thresfilter.h"

void calcDispls(int *src, int xsize, int ysize, int *displacements, int *sendCounts);

int main (int argc, char ** argv) {

  int error = 0;
  int* displacements;
  int* sendCounts;
  pixel* workData;
  MPI_Datatype pixelType;
  int workDataSize;

  /* MPI INIT*/
  int ierr = MPI_Init(&argc, &argv);
  int myId, numberProc;
  
  MPI_Comm_size(MPI_COMM_WORLD, &numberProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);
  
  int xsize, ysize, colmax;
  pixel src[MAX_PIXELS];
  double stime, etime;

  /* Take care of the arguments */

  if (argc != 3) {
    fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }

  /* read file */
  if(myId == 0)
  {
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
  	{
  	  error = 1;
  	}

    if (colmax > 255) {
      error = 1;
    }
  }

  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (error == 1) {
    fprintf(stderr, "Too large maximum color-component value\n");
    MPI_Finalize();
    exit(1);
  }

  //start time measure
  if(myId == 0)  
    stime = MPI_Wtime();

  /* Scatter the data */
  workDataSize = ysize * xsize / numberProc + 1;
  displacements = (int*) malloc(numberProc * sizeof(int));
  sendCounts = (int*) malloc(numberProc* sizeof(int));
  workData = (pixel*) malloc(workDataSize * sizeof(pixel));
  memset(workData, 0, workDataSize * sizeof(pixel));
  
  calcDispls(xsize, ysize, numberProc, displacements, sendCounts);
  constructPixelDataType(&pixelDataType);

  MPI_Scatterv(src, sendCounts, displacements, pixelType, workData, workDataSize, pixelType, 0, MPI_COMM_WORLD);

  if ( myId ==0 )
    printf("Has scattered the image, calling filter\n");

  thresfilter(workDataSize, workData, numberProc, myId);

  //  clock_gettime(CLOCK_REALTIME, &etime);

  //  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
  //	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

  MPI_Gatherv(workData, workDataSize, pixelType, src, sendCounts, displacements, pixelType, 0, MPI_COMM_WORLD);

  if(myId == 0){
    etime = MPI_Wtime();
    printf("Total time: %.2f", etime - stime);
  
    /* write result */
    printf("Writing output file\n");
    
    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      {
	MPI_Finalize();
	exit(1);
      }
  }
  MPI_Finalize();
  return(0);
}


void calcDispls(int xsize, int ysize, int nunProc, int *displacements, int *sendCounts)
{
  int currentDisplacement = 0;
  int sendCount;
  pixels = xsize * ysize / numProc;
  restPixels = xsize * ysize % numProc;
  
  for (int i = 0; i < numProc; i++) 
    {
      displacements[i] = currentDisplacement;
      
      sendCount = restPixels > 0 ? pixels + 1 : pixels;

      sendCounts[i] = sendCount;

      currentDisplacement += sendCount * sizeof(pixel);

      restPixels--;
    }
}

void constructPixelDataType(MPI_Datatype* pixelDataType) 
{
  pixel item;

  int block_lengths [] = {1, 1, 1};
  MPI_Datatype block_types [] = {MPI_CHAR, MPI_CHAR, MPI_CHAR};
  MPI_Aint start, displacement[3];

  MPI_Address(&item, &start);
  MPI_Address(&item.r, &displacement[0]);
  MPI_Address(&item.g, &displacement[1]);
  MPI_Address(&item.b, &displacement[2]);

  displacement[0] -= start;
  displacement[1] -= start;
  displacement[2] -= start;

  MPI_Type_struct(3, block_lengths, displacement, block_types, pixelDataType);

  MPI_Type_commit(pixelDataType);
}
