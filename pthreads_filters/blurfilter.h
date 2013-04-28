/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

struct thread_work_data{
    pixel *src, *target;
    int ystart, ystop;
    int radius;
};

struct thread_shared_data{
    int xsize, ysize;
    double *w;
};

struct thread_data{
    struct thread_shared_data* sharedData;
    struct thread_work_data* workData;
};



void* blurfilter(void *tParams);

#endif
