#!/bin/sh
#
#SBATCH -N 1 --tasks-per-node=2
#SBATCH -t 0:05:00
#

export OMP_NUM_THREADS=2

./ours
