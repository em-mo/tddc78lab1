#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "ppmio.h"
#include "thresfilter.h"

void calcDispls(int xsize, int ysize, int numThreads, int *displacements, int *writeCounts);

int main (int argc, char ** argv) {
    pthread_t *threads;
    int xsize, ysize, colmax, numThreads;
    pixel src[MAX_PIXELS];
    pixel target[MAX_PIXELS];
    struct timespec stime, etime;

    struct thread_shared_data shared_data;
    struct thread_data* thread_argument_data;

    int* displacements;
    int* sendCounts;

    /* Take care of the arguments */

    if (argc != 4) {
    	fprintf(stderr, "Usage: %s infile outfile numThreads\n", argv[0]);
    	exit(1);
    }

    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
    	fprintf(stderr, "Too large maximum color-component value\n");
    	exit(1);
    }

    numThreads  = atoi(argv[3]);

    threads = (pthread_t*) malloc(numThreads * sizeof(pthread_t));
    displacements = (int*) malloc(numThreads * sizeof(int));
    sendCounts = (int*) malloc(numThreads * sizeof(int));
    thread_argument_data = (struct thread_data*) malloc(numThreads * sizeof(struct thread_data));
  
    calcDispls(xsize, ysize, numThreads, displacements, sendCounts);

    shared_data.totalPixels = xsize*ysize;
    shared_data.sum = 0;
    pthread_mutex_init(&(shared_data.sumMutex), NULL);
    pthread_cond_init(&(shared_data.countThreadsCond), NULL);
    shared_data.sumThreadCount = numThreads;

    int t;
    for (t = 0; t < numThreads; t++)
    {
        thread_argument_data[t].shared_data = &shared_data;
        thread_argument_data[t].workData = &src[displacements[t]];
        thread_argument_data[t].workDataSize = sendCounts[t];

        if (pthread_create(&threads[t], NULL, thresfilter, (void*)&thread_argument_data[t])) 
        {
            printf("Error creating threads");
            for (;t >= 0; t--)
                pthread_cancel(threads[t]);
            exit(-1);
        }
    }

    for (t = 0; t < numThreads; t++) 
    {
        pthread_join(threads[t], NULL);
    }


    /* write result */
    printf("Writing output file\n");
    
    
    free(threads);
    free(displacements);
    free(sendCounts);

    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0) 
    {
        free(thread_argument_data);
        exit(1);
    }

    free(thread_argument_data);

    return(0);
}

void calcDispls(int xsize, int ysize, int numThreads, int *displacements, int *writeCounts)
{
    int currentDisplacement = 0;
    int writeCount, i, pixels, restPixels;
    pixels = xsize * ysize / numThreads;
    restPixels = (xsize * ysize) % numThreads;

    for (i = 0; i < numThreads; i++) 
    {
        displacements[i] = currentDisplacement;
        
        writeCount = restPixels > 0 ? pixels + 1 : pixels;
  
        writeCounts[i] = writeCount;
  
        currentDisplacement += writeCount;
  
        restPixels--;
    }
}
