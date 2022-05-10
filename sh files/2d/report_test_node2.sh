#!/bin/bash
#SBATCH -N 2
#SBATCH -C knl
#SBATCH -q debug
#SBATCH -J mpitest
#SBATCH -t 00:30:00

#modules
module load openmpi

# run the application
# fixed board size, freq, change total_rank
for sizex in 1000 2000
do
    for node in 2
    do
        for rank in 2 8 18 32 50 
        do
            echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi2d"
            srun -N $node --ntasks-per-node=$rank ./mpi2d -x $sizex -y 1000 -update 1 -t 128 -s 10
        done
    done
done