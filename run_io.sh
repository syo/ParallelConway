#!/bin/sh
#SBATCH -t 15
#SBATCH -N 8
#SBATCH -p debug
#SBATCH -D /gpfs/u/home/PCP6/PCP6fccj/scratch

srun -n 512 -N 8 --overcommit -o 1F_512R_1M.log ./assignment4 1 1048576
srun -n 512 -N 8 --overcommit -o 1F_512R_2M.log ./assignment4 1 2097152
srun -n 512 -N 8 --overcommit -o 1F_512R_4M.log ./assignment4 1 4194304
srun -n 512 -N 8 --overcommit -o 1F_512R_8M.log ./assignment4 1 8388608
srun -n 512 -N 8 --overcommit -o 1F_512R_16M.log ./assignment4 1 16777216
