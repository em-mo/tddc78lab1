#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=1
#SBATCH -t 0:00:10
#

mpprun ../parallel_filters/thresc ../images/im3.ppm out.ppm