#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=16
#SBATCH -t 0:00:10
#

mpprun ../parallel_filters/thresc ../images/im4.ppm out.ppm
