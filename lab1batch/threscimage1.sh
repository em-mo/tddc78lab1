#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:00:10
#

mpprun ../parallel_filters/thresc ../images/im1.ppm out.ppm
