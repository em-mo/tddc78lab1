#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "thresfilter.h"

void calcDispls(int xsize, int ysize, int numProc, int *displacements, int *sendCounts);
void constructPixelDataType(MPI_Datatype* pixelDataType);

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
  MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
  sendCounts = (int*) malloc(numberProc * sizeof(int));
  workData = (pixel*) malloc(workDataSize * sizeof(pixel));
  memset(workData, 0, workDataSize * sizeof(pixel));
  
  calcDispls(xsize, ysize, numberProc, displacements, sendCounts);
  MPI_Type_contiguous(3, MPI_CHAR, &pixelType);
  MPI_Type_commit(&pixelType);


  MPI_Scatterv(src, sendCounts, displacements, pixelType, workData, sendCounts[myId], pixelType, 0, MPI_COMM_WORLD);

  thresfilter(sendCounts[myId], workData, numberProc, myId, xsize * ysize);

  MPI_Gatherv(workData, sendCounts[myId], pixelType, src, sendCounts, displacements, pixelType, 0, MPI_COMM_WORLD);

  if(myId == 0){
    etime = MPI_Wtime();
    printf("Total time: %.6f", etime - stime);
  
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


void calcDispls(int xsize, int ysize, int numProc, int *displacements, int *sendCounts)
{
  int currentDisplacement = 0;
  int sendCount, i, pixels, restPixels;
  pixels = xsize * ysize / numProc;
  restPixels = (xsize * ysize) % numProc;
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);  

  for (i = 0; i < numProc; i++) 
    {
      displacements[i] = currentDisplacement;
      
      sendCount = restPixels > 0 ? pixels + 1 : pixels;

      sendCounts[i] = sendCount;

      currentDisplacement += sendCount;

      restPixels--;
      
      if(myId == 0){
	printf("displacement%d: %d\n", i, currentDisplacement);
      	printf("sendCount%d: %d\n", i, sendCount);
      }
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
