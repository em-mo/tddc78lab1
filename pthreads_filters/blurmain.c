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
pixel src[MAX_PIXELS];
pixel target[MAX_PIXELS];
double w[MAX_RAD];
int main (int argc, char ** argv) {

    int numThreads;
    int radius;
    
    int xsize, ysize, colmax;
    
    struct timespec stime, etime;
    
    int* displacements;
    int* writeCounts;

    

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

    printf("Has read the image, generating coefficients\n");

    /* filter */
    get_gauss_weights(radius, w);

    printf("Calling filter\n");

   // clock_gettime(CLOCK_REALTIME, &stime);

    
    printf("before calcDispls");
    calcDispls(xsize, ysize, numThreads, displacements, writeCounts);

    pthread_t threads[numThreads]; 
    struct thread_data thread_data_array[numThreads];
    struct thread_shared_data shared_data;
    shared_data.xsize = xsize;
    shared_data.ysize = ysize;
    shared_data.w = w;

    int t, ret;

    for(t=0; t < numThreads; t++)
    {   
        struct thread_work_data work_data;
        thread_data_array[t].workData = &work_data;
        thread_data_array[t].sharedData = &shared_data;

        work_data.ystart = displacements[t] / xsize;
        work_data.ystop = work_data.ystart + writeCounts[t] / xsize;
        work_data.src = src + displacements[t];
        work_data.target = target + displacements[t];
        work_data.radius = radius;
        ret = pthread_create(&threads[t], NULL, blurfilter, thread_data_array);
        if (ret) {
            printf("ERROR! Return code: %d\n", ret); 
            exit(-1); 
        }  
    }

    for(t=0; t < numThreads; t++)
    {
        pthread_join(threads[t], NULL);
    }

    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	    1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");
    
    if(write_ppm (argv[4], xsize, ysize, (char *)src) != 0)
      exit(1);


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
