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
  struct timespec stime, etime;

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
  	  MPI_Finalize();
  	  exit(1);
  	}

    if (colmax > 255) {
      error = 1;
    }
  }

  MPI_Bcast(&error, 1, MPI_Int, 0, MPI_COMM_WORLD);

  if (error == 1) {
    fprintf(stderr, "Too large maximum color-component value\n");
    MPI_Finalize();
    exit(1);
  }
  
  /* Scatter the data */
  workDataSize = ysize * xsize + 1;
  displacements = (int*) malloc(numberProc * sizeof(int));
  sendCounts = (int*) malloc(numberProc* sizeof(int));
  workData = (pixel*) malloc(workDataSize * sizeof(pixel));

  calcDispls(xsize, ysize, numberProc, displacements, sendCounts);
  constructPixelDataType(pixelDataType);

  MPI_Scatter(src, sendCounts, displacements, pixelType, workData, workDataSize, pixelType, 0, MPI_COMM_WORLD);

  printf("Has read the image, calling filter\n");

  clock_gettime(CLOCK_REALTIME, &stime);

  thresfilter(xsize, ysize, src);

  clock_gettime(CLOCK_REALTIME, &etime);

  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

  /* write result */
  printf("Writing output file\n");
    
  if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
    {
      MPI_Finalize();
      exit(1);
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

void constructPixelDataType(MPI_Datatype pixelDataType) 
{
  
}