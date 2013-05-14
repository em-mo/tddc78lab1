#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=1
#SBATCH -t 0:05:00
#

export OMP_NUM_THREADS=1

./ours
