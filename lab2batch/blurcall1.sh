#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=2
#SBATCH -t 0:00:10
#

../pthreads_filters/blurc 2 50 ../images/im1.ppm out.ppm
