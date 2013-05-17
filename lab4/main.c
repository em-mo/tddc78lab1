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

void doCommunication(vector<pcord_t> *neighbours, const int *neighbourRank, list<pcord_t *> *particles);
void initializeNeighbours(int neighbourCoords[][2], int *neighbourRanks, int *myCoords, int *dims);
void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds);
void initializeParticles(const int myId, const cord_t bounds, list<pcord_t *> *particles);
void changeLists(list<pcord_t *> *origin, const cord_t bounds, vector<pcord_t> *neighbours, int *neighbourRanks);
void resetLists(list<pcord_t *> *particles, list<pcord_t *> *collidedParticles, vector<pcord_t> *neighbours);
float randFloat();
void cleanUp(list<pcord_t *> *origin);

MPI_Datatype MPI_PCORD;
MPI_Comm gridComm;
int myId, numberProc;

int main(int argc, char **argv)
{
    //MPI
    int ierr = MPI_Init(&argc, &argv);
    int myCoords[2];

    myCoords[0] = 0;
    myCoords[1] = 0;


    //Cartesian Coordinates
    int dims[2];
    int periods[2];
    int reorder;

    cord_t bounds;

    int currentTimeStep = 0;

    cord_t wall;

    wall.x0 = 0;
    wall.y0 = 0;
    wall.x1 = BOX_HORIZ_SIZE;
    wall.y1 = BOX_VERT_SIZE;

    int neighbourCoords[4][2];
    int neighbourRanks[4];
    vector<pcord_t> neighbours[4];
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


    initializeNeighbours(neighbourCoords, neighbourRanks, myCoords, dims);
    initializeBounds(myCoords, dims, &bounds);
    initializeParticles(myId, bounds, &particles);

    if (myId == 0)
        printf("Initialization complete, beginning stepping\n");
    fflush(stdout);
    int t;
    float preassure = 0;
    while (currentTimeStep < MAX_STEPS)
    {
        int loopCount = 0;
        float tmpPreassure;
        for (list<pcord_t *>::iterator iter = particles.begin(); iter != particles.end(); ++iter)
        {
            pcord_t *p1 = *iter;
            tmpPreassure = 0;

            tmpPreassure = wall_collide(p1, wall);
            // if (myId == 0)
            //     printf("particle %d x %f y %f vx %f yx %f\n", ++loopCount, p1->x, p1->y, p1->vx, p1->vy);
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
        changeLists(&particles, bounds, neighbours, neighbourRanks);
        changeLists(&collidedParticles, bounds, neighbours, neighbourRanks);

        // COMMUNICATION
        doCommunication(neighbours, neighbourRanks, &particles);

        //RESET LISTS
        resetLists(&particles, &collidedParticles, neighbours);

        currentTimeStep += STEP_SIZE;
    }
    if (myId == 0)
        printf("Steps complete\n");

    resetLists(&particles, &collidedParticles, neighbours);

    printf("ID %d has %u particles\n", myId, (unsigned int)particles.size());

    //cleanUp(&particles);

    MPI_Type_free(&MPI_PCORD);
    MPI_Finalize();
}

void doCommunication(vector<pcord_t> *neighbours, const int *neighbourRank, list<pcord_t *> *particles)
{
    MPI_Request request;

    //for all neighbours
    for (int index = 0; index < 4; ++index)
    {
        if (neighbourRank[index] != -1)
        {
            int sendAmount = neighbours[index].size();
            MPI_Isend(&sendAmount, 1, MPI_INT, neighbourRank[index], TAG_AMOUNT, MPI_COMM_WORLD, &request);
            if (sendAmount != 0) 
            {
                // printf("Id %d, send %d to index %d\n", myId, sendAmount, neighbourRank[index]);
                MPI_Isend(&neighbours[index].front(), sendAmount, MPI_PCORD, neighbourRank[index], TAG_SEND, MPI_COMM_WORLD, &request);
            }
        }
    }

    //RECEIVE
    MPI_Status status;
    for (int index = 0; index < 4; ++index)
    {
        if (neighbourRank[index] != -1)
        {
            int recvAmount;
            MPI_Recv(&recvAmount, 1, MPI_INT, neighbourRank[index], TAG_AMOUNT, MPI_COMM_WORLD, &status);
            if (recvAmount != 0)
            {
                // printf("Id %d, recv %d to index %d\n", myId, recvAmount, neighbourRank[index]);
                pcord_t *recvBuffer = (pcord_t *)malloc(sizeof(pcord_t) * recvAmount);
                MPI_Recv(recvBuffer, recvAmount, MPI_PCORD, neighbourRank[index], TAG_SEND, MPI_COMM_WORLD, &status);
                for (int i = 0; i < recvAmount; i++)
                {
                    particles->push_back(&recvBuffer[i]);
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void resetLists(list<pcord_t *> *particles, list<pcord_t *> *collidedParticles, vector<pcord_t> *neighbours)
{
    while (!collidedParticles->empty())
    {
        particles->push_back(collidedParticles->back());
        collidedParticles->pop_back();
    }

    for (int index = 0; index < 4; ++index)
    {
        neighbours[index].clear();
    }
}

void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds)
{
    bounds->x0 = myCoords[0] * (BOX_HORIZ_SIZE / dims[0]);
    bounds->x1 = (myCoords[0] + 1) * (BOX_HORIZ_SIZE / dims[0]);
    bounds->y0 = myCoords[1] * (BOX_VERT_SIZE / dims[1]);
    bounds->y1 = (myCoords[1] + 1) * (BOX_VERT_SIZE / dims[1]);
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

        particle->vx = 2 * randFloat() * MAX_INITIAL_VELOCITY - MAX_INITIAL_VELOCITY;
        particle->vy = 2 * randFloat() * MAX_INITIAL_VELOCITY - MAX_INITIAL_VELOCITY;

        particles->push_back(particle);
    }
}

float randFloat()
{
    return (float)rand() / (float)RAND_MAX;
}

void changeLists(list<pcord_t *> *origin, const cord_t bounds, vector<pcord_t> *neighbours, int *neighbourRanks)
{
    pcord_t *particle;
    for (list<pcord_t *>::iterator iter = origin->begin(); iter != origin->end();)
    {
        particle = *iter;
        if (particle->x < bounds.x0 && neighbourRanks[INDEX_LEFT] != -1)
        {
            neighbours[INDEX_LEFT].push_back(*particle);
            //free(particle);
            iter = origin->erase(iter);
        }
        else if (particle->x > bounds.x1 && neighbourRanks[INDEX_RIGHT] != -1)
        {
            neighbours[INDEX_RIGHT].push_back(*particle);
            //free(particle);
            iter = origin->erase(iter);
        }
        else if (particle->y < bounds.y0 && neighbourRanks[INDEX_UP] != -1)
        {
            neighbours[INDEX_UP].push_back(*particle);
            //free(particle);
            iter = origin->erase(iter);
        }
        else if (particle->y > bounds.y1 && neighbourRanks[INDEX_DOWN] != -1)
        {
            neighbours[INDEX_DOWN].push_back(*particle);
            //free(particle);
            iter = origin->erase(iter);
        }
        else
        {
            iter++;
        }
    }
}

void initializeNeighbours(int neighbourCoords[][2], int *neighbourRanks, int *myCoords, int *dims)
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

void cleanUp(list<pcord_t *> *particles)
{
    for (list<pcord_t *>::iterator iter = particles->begin(); iter != particles->end(); ++iter)
    {
        free(*iter);
    }
}