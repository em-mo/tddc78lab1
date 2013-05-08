#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 0:00:10
#

../pthreads_filters/thresc ../images/im3.ppm out.ppm 4
