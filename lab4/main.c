#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <vector>
#include "definitions.h"

using namespace std;

#define NEIGHBOR_TOP = 0;
#define NEIGHBOR_RIGHT = 1;
#define NEIGHBOR_BOTTOM = 2;
#define NEIGHBOR_LEFT = 3;

#define MAX_STEPS 100

void initializeBounds(int *myCoords, int *dims, cord_t* bounds);
void initializeWalls(int *myCoords, int *dims, vector<cord_t*> *walls);

int main(int argc, char** argv)
{
  //MPI
  int ierr = MPI_Init(&argc, &argv);
  int myId, numberProc;
  int myCoords[2];

  //Cartesian Coordinates
  int dims[2];
  int periods[2];
  int reorder;
  MPI_Comm gridComm;
  
  cord_t bounds;

  int currentTimeStep = 0;
  
  vector<cord_t*> walls;
  vector<part_cord*> neighbourVectors[4];

 
  dims[0] = 0;
  dims[1] = 0;

  periods[0] = false;
  periods[1] = false;

  reorder = false;

  //Initialize MPI
  MPI_Comm_size(MPI_COMM_WORLD, &numberProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

  //Initialize Cartesian Coordinates
  MPI_Dims_create(numberProc, 2, dims);
  printf("x: %d, y: %d\n", dims[0], dims[1]);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &gridComm);

  MPI_Cart_rank(gridComm, myCoords, &myId);

  initializeBounds(myCoords, dims, &bounds);
  initializeWalls(myCoords, dims, &walls);
  
    
  MPI_Finalize();
}


void initializeBounds(int *myCoords, int *dims, cord_t* bounds)
{
  bounds->x0 = myCoords[0] * (BOX_HORIZ_SIZE / dims[0]);
  bounds->x1 = (myCoords[0]+1) * (BOX_HORIZ_SIZE / dims[0]);    
  bounds->y0 = myCoords[1] * (BOX_VERT_SIZE / dims[1]);
  bounds->y1 = (myCoords[1]+1) * (BOX_VERT_SIZE / dims[1]);
}

void initializeWalls(int *myCoords, int *dims, vector<cord_t*> *walls)
{
  cord_t* wall;
  if(myCoords[0] == 0)
    {
      wall = (cord_t*) malloc(sizeof(cord_t));
      wall->x0 = 0;
      wall->y0 = 0;
      wall->x1 = 0;
      wall->y1 = BOX_VERT_SIZE;
      walls->push_back(wall);
    }

  if(myCoords[1] == 0)
    {
      wall = (cord_t*) malloc(sizeof(struct cord));
      wall->x0 = 0;
      wall->y0 = 0;
      wall->x1 = BOX_HORIZ_SIZE;
      wall->y1 = 0;
      walls->push_back(wall);
    }  
  
  if(myCoords[0] == dims[0] - 1)
    {
      wall = (cord_t*) malloc(sizeof(struct cord));
      wall->x0 = BOX_HORIZ_SIZE;
      wall->y0 = 0;
      wall->x1 = BOX_HORIZ_SIZE;
      wall->y1 = BOX_VERT_SIZE;
      walls->push_back(wall);
    }
  
  if(myCoords[1] == dims[1] - 1)
    {
      wall = (cord_t*) malloc(sizeof(struct cord));
      wall->x0 = 0;
      wall->y0 = BOX_VERT_SIZE;
      wall->x1 = BOX_HORIZ_SIZE;
      wall->y1 = BOX_VERT_SIZE;
      walls->push_back(wall);
    }
}
