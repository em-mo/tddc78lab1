/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

#include <pthread.h>

/* NOTE: This structure must not be padded! */

typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

struct thread_shared_data {
    int totalPixels;
    unsigned int sum;
    pthread_mutex_t sumMutex;
    pthread_cond_t countThreadsCond;
    int sumThreadCount;
};

struct thread_data {
    struct thread_shared_data* shared_data;
    pixel* workData;
    int workDataSize;
};

void *thresfilter(void* params);

#endif
