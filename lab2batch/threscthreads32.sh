#!/bin/sh
#
#SBATCH -N 2 --tasks-per-node=16
#SBATCH -t 0:00:10
#

../pthreads_filters/thresc ../images/im3.ppm out.ppm 32
