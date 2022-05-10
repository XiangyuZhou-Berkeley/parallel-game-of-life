#!/bin/bash
#SBATCH -N 2
#SBATCH -C knl
#SBATCH -q debug
#SBATCH -J mpitest
#SBATCH -t 00:30:00

# # modules
# module load cmake
# module swap PrgEnv-intel PrgEnv-gnu
module load openmpi

# run the application
# fixed board size, freq, change total_rank 
for sizex in 1000 2000
do
    for node in 1
    do
        for rank in 1 2 4 8 9 16 25 30 32 36 49 60 64
        do
            echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi1d"
            srun -N $node --ntasks-per-node=$rank ./mpi1d -x $sizex -y 1000 -update 1 -t 128 -s 10
        done
    done
done
