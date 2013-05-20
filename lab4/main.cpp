#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

#define MAX_STEPS 100

void doCommunication(vector<pcord_t> *neighbourTransferData, const int *neighbourRank, list<pcord_t> *particles);
void initializeNeighbours(int neighbourCoords[][2], int *neighbourRanks, int *myCoords, int *dims);
void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds);
double initializeParticles(const int myId, const cord_t bounds, list<pcord_t> *particles);
void changeLists(list<pcord_t> *origin, const cord_t bounds, vector<pcord_t> *neighbourTransferData, int *neighbourRanks);
void resetLists(list<pcord_t> *particles, list<pcord_t> *collidedParticles, vector<pcord_t> *neighbourTransferData);
double randdouble();
void cleanUp(list<pcord_t> *particles);
void printParticles(list<pcord_t> *particles);
inline double calculatePressure(double pressure);
inline void safeIterDecrement(list<pcord_t>::iterator *iter, list<pcord_t> *iterList);


MPI_Datatype MPI_PCORD;
MPI_Comm gridComm;
int myId, numberProc;
int numberOfParticles, maxInitialVelocity, wallSideLengthX, wallSideLengthY;
double startTime, endTime;

int main(int argc, char **argv)
{
    //MPI
    int ierr = MPI_Init(&argc, &argv);
    int myCoords[2];

    myCoords[0] = 0;
    myCoords[1] = 0;

    numberOfParticles = INIT_NO_PARTICLES;
    maxInitialVelocity = MAX_INITIAL_VELOCITY;
    wallSideLengthX = BOX_HORIZ_SIZE;
    wallSideLengthY = BOX_VERT_SIZE;

    switch(argc)
    {
        case 4:
            wallSideLengthX = wallSideLengthY = atoi(argv[3]);
        case 3:
            maxInitialVelocity = atoi(argv[2]);
        case 2:
            numberOfParticles = atoi(argv[1]);
        default:
            break;
    }

    //Cartesian Coordinates
    int dims[2];
    int periods[2];
    int reorder;

    cord_t bounds;

    int currentTimeStep = 0;

    cord_t wall;

    wall.x0 = 0;
    wall.y0 = 0;
    wall.x1 = wallSideLengthX;
    wall.y1 = wallSideLengthY;

    int neighbourCoords[4][2];
    int neighbourRanks[4];
    vector<pcord_t> neighbourTransferData[4];
    list<pcord_t> particles;
    list<pcord_t> collidedParticles;

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

    MPI_Type_contiguous(PCORD_SIZE, MPI_DOUBLE, &MPI_PCORD);
    MPI_Type_commit(&MPI_PCORD);

    // if (myId == 0)
    //     printf("x: %d, y: %d\n", dims[0], dims[1]);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &gridComm);
    MPI_Cart_coords(gridComm, myId, 2, myCoords);

    double meanVelocity;

    startTime = MPI_Wtime();

    initializeNeighbours(neighbourCoords, neighbourRanks, myCoords, dims);
    initializeBounds(myCoords, dims, &bounds);
    meanVelocity = initializeParticles(myId, bounds, &particles);
    
    if(myId == 0)
    {
        printf("Dims %d %d\n", dims[0], dims[1]);
        printf("Particles %d  Max initial velocity %d  Wall length %d\n", numberOfParticles, maxInitialVelocity, wallSideLengthX);
        printf("Initialization complete, beginning stepping\n");
    }

    double t;
    double pressure = 0;
    double tmpPressure;
    bool hasCollided = false;
    while (currentTimeStep < MAX_STEPS)
    {
        for (list<pcord_t>::iterator iter = particles.begin(); iter != particles.end();)
        {
            pcord_t *p1 = &(*iter);
            tmpPressure = 0;

            tmpPressure = wall_collide(p1, wall);

            // Check for particle collisions
            if (tmpPressure == 0)
            {
                list<pcord_t>::iterator iter2(iter);
                for (++iter2; iter2 != particles.end(); ++iter2)
                {
                    pcord_t *p2 = &(*iter2);

                    t = collide(p1, p2);
                    if (t != -1)
                    {
                        interact(p1, p2, t);
                        collidedParticles.push_back(*iter);
                        collidedParticles.push_back(*iter2);
                        particles.erase(iter2);
                        iter = particles.erase(iter);
                        hasCollided = true;
                        break;
                    }
                }
            }
            else //collided with wall
            {
                pressure += tmpPressure;
                collidedParticles.push_back(*iter);
                iter = particles.erase(iter);
                hasCollided = true;
            }

            if (hasCollided)
                hasCollided = false;
            else
                ++iter;
            

        }

        for (list<pcord_t>::iterator it = particles.begin(); it != particles.end(); ++it)
        {
            feuler(&(*it), STEP_SIZE);
        }
        //Move particles to neighbour lists that should be moved to new porocess.
        changeLists(&particles, bounds, neighbourTransferData, neighbourRanks);
        changeLists(&collidedParticles, bounds, neighbourTransferData, neighbourRanks);

        // COMMUNICATION
        doCommunication(neighbourTransferData, neighbourRanks, &particles);

        //RESET LISTS
        resetLists(&particles, &collidedParticles, neighbourTransferData);

        currentTimeStep += STEP_SIZE;
    }

    if (myId == 0)
        printf("Steps complete\n");

    double totalPressure;
    MPI_Reduce(&pressure, &totalPressure, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    endTime = MPI_Wtime();

    int area;

    double rt;

    if (myId == 0)
    {
        pressure = calculatePressure(totalPressure);
        rt = pressure * (wallSideLengthX * wallSideLengthY) / (numberOfParticles * numberProc);
        printf("Wall pressure is %f with RT %f in %f seconds\n", calculatePressure(totalPressure), rt, endTime - startTime);
    }
    MPI_Type_free(&MPI_PCORD);
    MPI_Finalize();
}

void doCommunication(vector<pcord_t> *neighbourTransferData, const int *neighbourRank, list<pcord_t> *particles)
{
    MPI_Request request;

    //for all neighbour SEND
    for (int index = 0; index < 4; ++index)
    {
        if (neighbourRank[index] != -1)
        {
            int sendAmount = neighbourTransferData[index].size();
            MPI_Isend(&sendAmount, 1, MPI_INT, neighbourRank[index], TAG_AMOUNT, MPI_COMM_WORLD, &request);
            if (sendAmount != 0) 
            {
                // printf("Id %d, send %d to index %d\n", myId, sendAmount, neighbourRank[index]);
                MPI_Isend(&(neighbourTransferData[index].front()), sendAmount, MPI_PCORD, neighbourRank[index], TAG_SEND, MPI_COMM_WORLD, &request);
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
                    particles->push_back(recvBuffer[i]);
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void resetLists(list<pcord_t> *particles, list<pcord_t> *collidedParticles, vector<pcord_t> *neighbourTransferData)
{
    while (!collidedParticles->empty())
    {
        particles->push_back(collidedParticles->back());
        collidedParticles->pop_back();
    }

    for (int index = 0; index < 4; ++index)
    {
        neighbourTransferData[index].clear();
    }
}

void initializeBounds(const int *myCoords, const int *dims, cord_t *bounds)
{
    bounds->x0 = myCoords[0] * (wallSideLengthX / dims[0]);
    bounds->x1 = (myCoords[0] + 1) * (wallSideLengthX / dims[0]);
    bounds->y0 = myCoords[1] * (wallSideLengthY / dims[1]);
    bounds->y1 = (myCoords[1] + 1) * (wallSideLengthY / dims[1]);
}

double initializeParticles(const int myId, const cord_t bounds, list<pcord_t> *particles)
{
    double meanVelocityX = 0;
    double meanVelocityY = 0;

    int amount = numberOfParticles + 1;
    pcord_t particle;
    srand(time(NULL) * (myId + 1));

    int boundsWidth = bounds.x1 - bounds.x0;
    int boundsHeight = bounds.y1 - bounds.y0;
	double r;
	double angle;

    while (--amount)
    {
        particle.x = randdouble() * boundsWidth + bounds.x0;
        particle.y = randdouble() * boundsHeight + bounds.y0;


        r = randdouble() * maxInitialVelocity;
		angle = randdouble() * 2 * PI;
		particle.vx = r * cos(angle);
		particle.vy = r * sin(angle);

        meanVelocityX += fabs(particle.vx);
        meanVelocityY += fabs(particle.vy);

        particles->push_back(particle);
	}

	return sqrtf(meanVelocityX * meanVelocityX + meanVelocityY * meanVelocityY) / numberOfParticles;
}

double randdouble()
{
    return (double)rand() / (double)RAND_MAX;
}

void changeLists(list<pcord_t> *origin, const cord_t bounds, vector<pcord_t> *neighbourTransferData, int *neighbourRanks)
{
    pcord_t particle;
    for (list<pcord_t>::iterator iter = origin->begin(); iter != origin->end();)
    {
        particle = *iter;
        if (particle.x < bounds.x0 && neighbourRanks[INDEX_LEFT] != -1)
        {
            neighbourTransferData[INDEX_LEFT].push_back(particle);
            // free(particle);
            iter = origin->erase(iter);
        }
        else if (particle.x > bounds.x1 && neighbourRanks[INDEX_RIGHT] != -1)
        {
            neighbourTransferData[INDEX_RIGHT].push_back(particle);
            // free(particle);
            iter = origin->erase(iter);
        }
        else if (particle.y < bounds.y0 && neighbourRanks[INDEX_UP] != -1)
        {
            neighbourTransferData[INDEX_UP].push_back(particle);
            // free(particle);
            iter = origin->erase(iter);
        }
        else if (particle.y > bounds.y1 && neighbourRanks[INDEX_DOWN] != -1)
        {
            neighbourTransferData[INDEX_DOWN].push_back(particle);
            // free(particle);
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

void cleanUp(list<pcord_t> *particles)
{
    for (list<pcord_t>::iterator iter = particles->begin(); iter != particles->end(); ++iter)
    {
        // free(*iter);
    }
}

void printParticles(list<pcord_t> *particles) 
{
    for (list<pcord_t>::iterator iter = particles->begin(); iter != particles->end(); ++iter)
    {
        printf("ID %d particle x %f y %f vx %f vy %f\n", myId, iter->x, iter->y, iter->vx, iter->vy);
    }
}

inline double calculatePressure(double pressure)
{
    return pressure / ((wallSideLengthX * 2 + wallSideLengthY * 2) * MAX_STEPS);
}

inline void safeIterDecrement(list<pcord_t>::iterator *iter, list<pcord_t> *iterList)
{
    if (*iter != iterList->begin())
        --(*iter);
}