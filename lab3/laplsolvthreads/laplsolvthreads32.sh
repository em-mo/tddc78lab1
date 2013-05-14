#!/bin/sh
#
#SBATCH -N 2 --tasks-per-node=16
#SBATCH -t 0:05:00
#

export OMP_NUM_THREADS=32

./ours
