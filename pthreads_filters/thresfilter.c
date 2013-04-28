#include "thresfilter.h"

void *thresfilter(void *params)
{
#define uint unsigned int 

    uint sum, i, average, psum;
    
    struct thread_data* thread_argument_data = (struct thread_data*) params;

    pixel* src = thread_argument_data->workData;
    const int workDataSize = thread_argument_data->workDataSize;
    unsigned int* finalSum = &(thread_argument_data->shared_data->sum);
    pthread_mutex_t* sumMutex = &(thread_argument_data->shared_data->sumMutex);
    const int totalPixels = thread_argument_data->shared_data->totalPixels;
    pthread_cond_t* countThreadsCond = &(thread_argument_data->shared_data->countThreadsCond);
    int* sumThreadCount = &(thread_argument_data->shared_data->sumThreadCount);

    for(i = 0, sum = 0; i < workDataSize; i++) 
    {
        sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    }
  
    pthread_mutex_lock(sumMutex);
    *finalSum += sum;
    *sumThreadCount -= 1;
    if (*sumThreadCount == 0)
      pthread_cond_broadcast(countThreadsCond);
    while (*sumThreadCount != 0)
      pthread_cond_wait(countThreadsCond, sumMutex);

    pthread_mutex_unlock(sumMutex);

    average = *finalSum / totalPixels;
  
    for(i = 0; i < workDataSize; i++) 
    {
        psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
        if(average > psum) 
        {
            src[i].r = src[i].g = src[i].b = 0;
        }
        else 
        {
            src[i].r = src[i].g = src[i].b = 255;
        }
    }
}
