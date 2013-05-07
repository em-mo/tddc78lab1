#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=1
#SBATCH -t 0:00:10
#

../pthreads_filters/blurc 1 50 ../images/im3.ppm out.ppm
