#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <list>
#include <vector>
#include "definitions.h"

using namespace std;

#define INDEX_UP 0
#define INDEX_RIGHT 1
#define INDEX_DOWN 2
#define INDEX_LEFT 3

#define MAX_STEPS 100


void initializeBounds(const int *myCoords, const int *dims, cord_t* bounds);
void initializeWalls(const int *myCoords, const int *dims, vector<cord_t*> *walls);
void  initializeParticles(const int myId, const cord_t bounds, list<pcord_t*> *particles);
void changeLists(list<pcord_t*> l, const cord_t bounds, list<pcord_t*> *neighbours);
float randFloat();

int main(int argc, char** argv)
{
  //MPI
  int ierr = MPI_Init(&argc, &argv);
  int myId, numberProc;
  int myCoords[2];

  myCoords[0] = 0;
  myCoords[1] = 0;


  //Cartesian Coordinates
  int dims[2];
  int periods[2];
  int reorder;
  MPI_Comm gridComm;
  
  cord_t bounds;

  int currentTimeStep = 0;
  
  vector<cord_t*> walls;
  list<pcord_t*> neighbours[4];
  list<pcord_t*> particles;
  list<pcord_t*> collidedParticles;
 
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

  if (myId == 0)
    printf("x: %d, y: %d\n", dims[0], dims[1]);

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &gridComm);
  MPI_Cart_coords(gridComm, myId, 2, myCoords);

  initializeBounds(myCoords, dims, &bounds);
  initializeWalls(myCoords, dims, &walls);
  initializeParticles(myId, bounds, &particles);

  float preassure = 0;
  while(currentTimeStep < MAX_STEPS)
    {
      int t;
      float tmpPreassure;
      for(list<pcord_t*>::iterator iter = particles.begin(); iter != particles.end(); iter++)
	{
	  pcord_t* p1 = *iter;
	  tmpPreassure = 0;

	  //Check collision with walls
	  for(vector<cord_t*>::iterator wIter = walls.begin(); wIter != walls.end() && tmpPreassure == 0; wIter++)
	    {
	      tmpPreassure = wall_collide(p1, **wIter);
	    }

	  //Check for particle collisions
	  if(tmpPreassure == 0)
	    {
        list<pcord_t*>::iterator iter2 = iter;
	      for(iter2++; iter2 != particles.end(); iter2++)
		{
		  pcord_t* p2 = *iter2;

		  t = collide(p1,p2);
		  if(t != -1)
		    {
		      interact(p1, p2, t);
		      collidedParticles.push_back(*iter);
		      collidedParticles.push_back(*iter2);
		      iter = particles.erase(iter);
		      particles.erase(iter2);
		      --iter;
		      break;
		    }
		}
	    }
	  else //collided with wall
	    {
	      preassure += tmpPreassure;
	      collidedParticles.push_back(*iter);
	      iter = particles.erase(iter);
	      --iter;
	    }
	}
      
      //Move particles to neighbour lists that should be moved to new porocess.
      changeLists(particles, bounds, neighbours);
      changeLists(collidedParticles, bounds, neighbours);
      
      // COMMUNICATION
      

      //RESET LISTS

      currentTimeStep += STEP_SIZE;
    }

  MPI_Finalize();
}


void initializeBounds(const int *myCoords, const int *dims, cord_t* bounds)
{
  bounds->x0 = myCoords[0] * (BOX_HORIZ_SIZE / dims[0]);
  bounds->x1 = (myCoords[0]+1) * (BOX_HORIZ_SIZE / dims[0]);    
  bounds->y0 = myCoords[1] * (BOX_VERT_SIZE / dims[1]);
  bounds->y1 = (myCoords[1]+1) * (BOX_VERT_SIZE / dims[1]);
}

void initializeWalls(const int *myCoords, const int *dims, vector<cord_t*> *walls)
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

void  initializeParticles(const int myId, const cord_t bounds, list<pcord_t*> *particles)
{
  int amount = INIT_NO_PARTICLES + 1;
  pcord_t* particle;
  srand(time(NULL) * (myId + 1));

  int boundsWidth = bounds.x1 - bounds.x0;
  int boundsHeight = bounds.y1 - bounds.y0;
  while(--amount)
    {
      particle = (pcord_t*) malloc(sizeof(pcord_t));
      particle->x = randFloat() * boundsWidth + bounds.x0;
      particle->y = randFloat() * boundsHeight + bounds.y0;

      particle->vx = randFloat() * MAX_INITIAL_VELOCITY;
      particle->vy = randFloat() * MAX_INITIAL_VELOCITY;
      
      particles->push_back(particle);
    }
}

float randFloat()
{
  return (float)rand()/(float)RAND_MAX;
}

void changeLists(list<pcord_t*> l, const cord_t bounds, list<pcord_t*> *neighbours)
{
  pcord_t* particle;
  for(list<pcord_t*>::iterator iter = l.begin(); iter != l.end();)
    {
      particle = *iter;
      if(particle->x < bounds.x0)
	{
	  neighbours[INDEX_LEFT].push_back(particle);
	  iter = l.erase(iter);
	}
      else if(particle->x > bounds.x1)
	{
	  neighbours[INDEX_RIGHT].push_back(particle);
	  iter = l.erase(iter);
	}
      else if(particle->y < bounds.y0)
	{
	  neighbours[INDEX_UP].push_back(particle);
	  iter = l.erase(iter);
	}
      else if(particle->y > bounds.y1)
	{
	  neighbours[INDEX_DOWN].push_back(particle);
	  iter = l.erase(iter);
	}
      else
	{
	  iter++;
	}
    }
}
