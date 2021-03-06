/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

#include <mpi.h>
/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

void thresfilter(const int nump, pixel* src, int numberProc, int myId, int totalPixels);

#endif
