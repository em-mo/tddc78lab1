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

#define TAG_AMOUNT 444
#define TAG_SEND 555

#define PCORD_SIZE 4

#define MAX_STEPS 10


void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds);
void initializeWalls(const int *myCoords, const int *dims, vector<cord_t *> *walls);
void  initializeParticles(const int myId, const cord_t bounds, list<pcord_t *> *particles);
void changeLists(list<pcord_t *> *origin, const cord_t bounds, list<pcord_t *> *neighbours);
void resetLists(list<pcord_t *> *particles, list<pcord_t *> *collidedParticles);
float randFloat();

int main(int argc, char **argv)
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
    MPI_datatype MPI_PCORD;
    
    cord_t bounds;

    int currentTimeStep = 0;

    cord_t wall;

    wall.x0 = 0;
    wall.y0 = 0;
    wall.x1 = BOX_HORIZ_SIZE;
    wall.y1 = BOX_VERT_SIZE;

    vector<cord_t *> walls;
    int neighbourCoords[4][2];
    int neighbourRanks[4];
    vector<pcord_t *> neighbours[4];
    list<pcord_t *> particles;
    list<pcord_t *> collidedParticles;

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

    MPI_Type_contiguous(PCORD_SIZE, MPI_FLOAT, &MPI_PCORD);
    MPI_Type_commit(&MPI_PCORD);

    // if (myId == 0)
    //     printf("x: %d, y: %d\n", dims[0], dims[1]);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &gridComm);
    MPI_Cart_coords(gridComm, myId, 2, myCoords);
    

    initializeNeighbours(neighbourCoords, myCoords);
    initializeBounds(myCoords, dims, &bounds);
    initializeWalls(myCoords, dims, &walls);
    initializeParticles(myId, bounds, &particles);

    printf("Id: %u  walls: %u\n", myId, (unsigned int)walls.size());

    int t;
    float preassure = 0;
    while (currentTimeStep < MAX_STEPS)
    {
        float tmpPreassure;
        for (list<pcord_t *>::iterator iter = particles.begin(); iter != particles.end(); ++iter)
        {
            pcord_t *p1 = *iter;
            tmpPreassure = 0;

            tmpPreassure = wall_collide(p1, wall);

            preassure += tmpPreassure;
            //Check for particle collisions
            if (tmpPreassure == 0)
            {
                list<pcord_t *>::iterator iter2 = iter;
                for (iter2++; iter2 != particles.end(); ++iter2)
                {
                    pcord_t *p2 = *iter2;

                    t = collide(p1, p2);
                    if (t != -1)
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

        for (list<pcord_t *>::iterator iter = particles.begin(); iter != particles.end(); ++iter) 
        {
            feuler(*iter, STEP_SIZE);
        }

        //Move particles to neighbour lists that should be moved to new porocess.
        changeLists(&particles, bounds, neighbours);
        changeLists(&collidedParticles, bounds, neighbours);

        // COMMUNICATION
        printf("ID: %u   collidedParticles: %u \n", myId, (unsigned int)collidedParticles.size());
        printf("ID: %u   other   Particles: %u \n", myId, (unsigned int)particles.size());

        //RESET LISTS
        resetLists(&particles, &collidedParticles);

        currentTimeStep += STEP_SIZE;
    }
    
    MPI_Type_free(&MPI_PCORD);
    MPI_Finalize();
}

void doCommunication(list<pcord_t *> *neighbours, int* neighbourRank)
{
    MPI_Request request;

    //for all neighbours
    for(int index = 0; index < 4; index++)
    {
	if (neighbourRank[index] != -1)
	{
	    int sendAmount = neighbours[index].size();
	    MPI_Isend(sendAmount, 1, MPI_INT, neighbourRank[index], TAG_AMOUNT, MPI_COMM_WORLD, &request);
	    if (sendAmount != 0)
		MPI_Isend(&neighbours[index].front(), sendAmount, MPI_PCORD, neighbourRank[index], TAG_SEND, MPI_COMM_WORLD, &request);
	}
    }

    //RECEIVE
    MPI_Status status;
    float* recvBuf = malloc(sizeof(float) * COMM_BUFFER_SIZE);
    for(int index = 0; index < 4; index++)
    {
	if (neighbourRank[index] != -1)
	{
	    int recvAmount;
	    MPI_Recv(&recvAmount, 1, MPI_INT, neighbourRank[index], TAG_AMOUNT, MPI_COMM_WORLD, &status);
	    if (recvAmount != 0)
	    {
		MPI_Recv(&neighbours[index].front(), sendAmount, MPI_PCORD, neighbourRank[index], TAG_SEND, MPI_COMM_WORLD, &request);
	    }
	}
    }    
}

void resetLists(list<pcord_t *> *particles, list<pcord_t *> *collidedParticles)
{
    while (!collidedParticles->empty())
    {
        particles->push_back(collidedParticles->back());
        collidedParticles->pop_back();
    }
}

void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds)
{
    bounds->x0 = myCoords[0] * (BOX_HORIZ_SIZE / dims[0]);
    bounds->x1 = (myCoords[0] + 1) * (BOX_HORIZ_SIZE / dims[0]);
    bounds->y0 = myCoords[1] * (BOX_VERT_SIZE / dims[1]);
    bounds->y1 = (myCoords[1] + 1) * (BOX_VERT_SIZE / dims[1]);
}

void initializeWalls(const int *myCoords, const int *dims, vector<cord_t *> *walls)
{
    cord_t *wall;
    if (myCoords[0] == 0)
    {
        wall = (cord_t *) malloc(sizeof(cord_t));
        wall->x0 = 0;
        wall->y0 = 0;
        wall->x1 = 0;
        wall->y1 = BOX_VERT_SIZE;
        walls->push_back(wall);
    }

    if (myCoords[1] == 0)
    {
        wall = (cord_t *) malloc(sizeof(struct cord));
        wall->x0 = 0;
        wall->y0 = 0;
        wall->x1 = BOX_HORIZ_SIZE;
        wall->y1 = 0;
        walls->push_back(wall);
    }

    if (myCoords[0] == dims[0] - 1)
    {
        wall = (cord_t *) malloc(sizeof(struct cord));
        wall->x0 = BOX_HORIZ_SIZE;
        wall->y0 = 0;
        wall->x1 = BOX_HORIZ_SIZE;
        wall->y1 = BOX_VERT_SIZE;
        walls->push_back(wall);
    }

    if (myCoords[1] == dims[1] - 1)
    {
        wall = (cord_t *) malloc(sizeof(struct cord));
        wall->x0 = 0;
        wall->y0 = BOX_VERT_SIZE;
        wall->x1 = BOX_HORIZ_SIZE;
        wall->y1 = BOX_VERT_SIZE;
        walls->push_back(wall);
    }
}

void  initializeParticles(const int myId, const cord_t bounds, list<pcord_t *> *particles)
{
    int amount = INIT_NO_PARTICLES + 1;
    pcord_t *particle;
    srand(time(NULL) * (myId + 1));

    int boundsWidth = bounds.x1 - bounds.x0;
    int boundsHeight = bounds.y1 - bounds.y0;
    while (--amount)
    {
        particle = (pcord_t *) malloc(sizeof(pcord_t));
        particle->x = randFloat() * boundsWidth + bounds.x0;
        particle->y = randFloat() * boundsHeight + bounds.y0;

        particle->vx = randFloat() * MAX_INITIAL_VELOCITY;
        particle->vy = randFloat() * MAX_INITIAL_VELOCITY;

        particles->push_back(particle);
    }
}

float randFloat()
{
    return (float)rand() / (float)RAND_MAX;
}

void changeLists(list<pcord_t *> *origin, const cord_t bounds, list<pcord_t *> *neighbours)
{
    pcord_t *particle;
    for (list<pcord_t *>::iterator iter = origin->begin(); iter != origin->end();)
    {
        particle = *iter;
        if (particle->x < bounds.x0)
        {
            neighbours[INDEX_LEFT].push_back(particle);
            iter = origin->erase(iter);
        }
        else if (particle->x > bounds.x1)
        {
            neighbours[INDEX_RIGHT].push_back(particle);
            iter = origin->erase(iter);
        }
        else if (particle->y < bounds.y0)
        {
            neighbours[INDEX_UP].push_back(particle);
            iter = origin->erase(iter);
        }
        else if (particle->y > bounds.y1)
        {
            neighbours[INDEX_DOWN].push_back(particle);
            iter = origin->erase(iter);
        }
        else
        {
            iter++;
        }
    }
}

void initializeNeighbours(MPI_Comm gridComm, int** neighbourCoords, int* neighbourRanks, int* myCoords, int* dims)
{
    for (int i = 0; i < 4; ++i) 
    {
	neighbourRanks[i] = -1;
    }

    neighbourCoords[INDEX_UP][0] = myCoords[0]; 
    neighbourCoords[INDEX_UP][1] = myCoords[1] - 1;
    if (neighbourCoords[INDEX_UP][1] != -1)
	MPI_Cart_rank(gridComm, neighbourCoords[INDEX_UP], &neighbourRanks[INDEX_UP]);
    
    neighbourCoords[INDEX_RIGHT][0] = myCoords[0] + 1; 
    neighbourCoords[INDEX_RIGHT][1] = myCoords[1];
    if (neighbourCoords[INDEX_RIGHT][0] != dims[0])
	MPI_Cart_rank(gridComm, neighbourCoords[INDEX_RIGHT], &neighbourRanks[INDEX_RIGHT]);

    neighbourCoords[INDEX_DOWN][0] = myCoords[0]; 
    neighbourCoords[INDEX_DOWN][1] = myCoords[1] + 1;
    if (neighbourCoords[INDEX_DOWN][1] != dims[1])
	MPI_Cart_rank(gridComm, neighbourCoords[INDEX_DOWN], &neighbourRanks[INDEX_DOWN]);

    neighbourCoords[INDEX_LEFT][0] = myCoords[0] - 1; 
    neighbourCoords[INDEX_LEFT][1] = myCoords[1];
    if (neighbourCoords[INDEX_LEFT][0] != -1)
	MPI_Cart_rank(gridComm, neighbourCoords[INDEX_LEFT], &neighbourRanks[INDEX_LEFT]);
}
