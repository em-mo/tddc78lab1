#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:01:30
#

../pthreads_filters/blurc 8 50 ../images/im4.ppm out.ppm
