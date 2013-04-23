#include "thresfilter.h"

void thresfilter(const int nump, pixel* src, int numberProc, int myId ){
#define uint unsigned int 

  uint sum, totalSum, i, psum, nump;

  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }
  
  MPI_Allreduce(&sum, &totalSum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  totalSum /= nump;
  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(totalSum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}

