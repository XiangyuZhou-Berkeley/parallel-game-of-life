#!/bin/bash
#SBATCH -N 2
#SBATCH -C knl
#SBATCH -q debug
#SBATCH -J mpitest
#SBATCH -t 00:30:00

#modules
module load openmpi

# run the application
# Weak scaling
for node in 1
do
    for rank in 1 2 4 8 9 16 25 32 36 40 50 60 64
    do
        sizex=$(( 1000*$((node * rank)) ))
        echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi1d"
        srun -N $node --ntasks-per-node=$rank ./mpi1d -x $sizex -y 1000 -update 1 -t 128 -s 10
    done
done

for node in 2
do
    for rank in 1 2 4 8 9 16 18 25 32 36 40 50 60 64
    do
        sizex=$(( 1000*$((node * rank)) ))
        echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi1d"
        srun -N $node --ntasks-per-node=$rank ./mpi1d -x $sizex -y 1000 -update 1 -t 128 -s 10
    done
done