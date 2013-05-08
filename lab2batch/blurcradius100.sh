#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:00:20
#

../pthreads_filters/blurc 8 100 ../images/im3.ppm out.ppm
