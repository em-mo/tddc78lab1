#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

void calcDispls(int xsize, int ysize, int numThreads, int *displacements, int *writeCounts);

#define MAX_RAD 1000

double w[MAX_RAD];

int main (int argc, char ** argv) {

    int numThreads;
    int radius;
    
    int xsize, ysize, colmax;
    
    struct timespec stime, etime;
    
    int* displacements;
    int* writeCounts;

    pixel* src = (pixel*) malloc(MAX_PIXELS * sizeof(pixel));
    pixel* target = (pixel*) malloc(MAX_PIXELS * sizeof(pixel));

    printf("Pthreads blur Threads: %s, radius: %s, file: %s\n", argv[1], argv[2], argv[3]);

    /* Take care of the arguments */

    if (argc != 5) {
        fprintf(stderr, "Usage: %s numthreads radius infile outfile\n", argv[0]);
	    exit(1);
    }

    numThreads = atoi(argv[1]);
    radius = atoi(argv[2]);

    if((radius > MAX_RAD) || (radius < 1)) {
	    fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
	    exit(1);
    }

    /* read file */
    if(read_ppm (argv[3], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
	fprintf(stderr, "Too large maximum color-component value\n");
	exit(1);
    }

    memcpy(target, src, MAX_PIXELS * sizeof(pixel));

    /* filter */
    get_gauss_weights(radius, w);


    clock_gettime(CLOCK_REALTIME, &stime);

    displacements = (int*) malloc(numThreads * sizeof(int));
    writeCounts = (int*) malloc(numThreads * sizeof(int));

    /* Calculate displacements */
    calcDispls(xsize, ysize, numThreads, displacements, writeCounts);

    /*Initiate shared data for threads*/
    pthread_t threads[numThreads]; 
    struct thread_data thread_data_array[numThreads];
    struct thread_shared_data shared_data;
    shared_data.xsize = xsize;
    shared_data.ysize = ysize;
    shared_data.w = w;

    int t, ret;

    /* Initiate work data for thread and create new threads */
    for(t=0; t < numThreads; t++)
    {   
        struct thread_work_data* work_data = (struct thread_work_data*) malloc(sizeof(struct thread_work_data));
        thread_data_array[t].workData = work_data;
        thread_data_array[t].sharedData = &shared_data;

        work_data->ystart = displacements[t] / xsize;
        work_data->ystop = work_data->ystart + writeCounts[t] / xsize;
        work_data->src = src;
        work_data->target = target;
        work_data->radius = radius;
        work_data->myId = t;

        ret = pthread_create(&threads[t], NULL, blurfilter, &thread_data_array[t]);

        if (ret) {
            printf("ERROR! Return code: %d\n", ret); 
            exit(-1); 
        }  
    }

    /* Join threads */
    for(t=0; t < numThreads; t++)
    {
        pthread_join(threads[t], NULL);
    }

    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	    1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    
    if(write_ppm (argv[4], xsize, ysize, (char *)target) != 0)
      exit(1);

    free(displacements);
    free(writeCounts);
    free(src);
    free(target);

    for(t=0; t < numThreads; t++)
    {
        free(thread_data_array[t].workData);
    }

    return(0);
}

void calcDispls(int xsize, int ysize, int numThreads, int *displacements, int *writeCounts)
{
    int currentDisplacement = 0;
    int writeCount, i, rows, restRows;
    rows = ysize / numThreads;
    restRows = ysize % numThreads;

    for (i = 0; i < numThreads; i++) 
    {
        displacements[i] = currentDisplacement;
        writeCount = restRows > 0 ? rows + 1 : rows;
        writeCount *= xsize;
        writeCounts[i] = writeCount;

        currentDisplacement += writeCount;

        restRows--;
    }
}
