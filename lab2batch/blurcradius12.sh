#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:00:10
#

../pthreads_filters/blurc 8 12 ../images/im3.ppm out.ppm
