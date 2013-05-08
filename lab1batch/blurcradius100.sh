#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=8
#SBATCH -t 0:00:20
#

mpprun ../parallel_filters/blurc 100 ../images/im3.ppm out.ppm
