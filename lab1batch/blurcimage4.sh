#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:01:30
#

mpprun ../parallel_filters/blurc 50 ../images/im4.ppm out.ppm
