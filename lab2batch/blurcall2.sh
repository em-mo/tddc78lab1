#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 0:00:10
#

../pthreads_filters/blurc 4 50 ../images/im2.ppm out.ppm
