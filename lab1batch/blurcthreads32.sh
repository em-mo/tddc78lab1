#!/bin/sh
#
#SBATCH -N 2 --tasks-per-node=16
#SBATCH -t 0:00:10
#

mpprun ../parallel_filters/blurc 50 ../images/im3.ppm out.ppm
