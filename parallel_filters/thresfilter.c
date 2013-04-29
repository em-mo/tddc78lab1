#include "thresfilter.h"

void thresfilter(const int nump, pixel* src, int numberProc, int myId, int totalPixels){
#define uint unsigned int 

  uint sum, i, psum;
  uint totalSum = 0;

  double avg;
  // Calculate a private sum, get global total sum and then take the average
  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }
  
  MPI_Allreduce(&sum, &totalSum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  avg = totalSum / totalPixels;

  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(avg > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}

