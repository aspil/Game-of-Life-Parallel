#!/bin/bash

# Max VM size #
#PBS -l pvmem=2G

# Max Wall time, Example 1 Minute #
#PBS -l walltime=00:01:00

# How many nodes and tasks per node, Example 2 nodes with 8 tasks each => 16 tasts #
#PBS -l select=1:ncpus=8:mpiprocs=4

# Only this job uses the chosen nodes
#PBS -l place=excl

# Which Queue to use, DO NOT CHANGE #
#PBS -q workq

# JobName #
#PBS -N myJob_4

#Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR
# Run executable #
mpirun -n 4 game.x
