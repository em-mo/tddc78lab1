#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=2
#SBATCH -t 0:00:10
#

../pthreads_filters/thresc ../images/im1.ppm out.ppm 2
