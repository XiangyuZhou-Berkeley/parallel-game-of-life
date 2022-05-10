#!/bin/bash
#SBATCH -N 2
#SBATCH -C knl
#SBATCH -q debug
#SBATCH -J mpitest
#SBATCH -t 00:30:00

#modules
module load openmpi
# run the application
# fixed board size, change freq, total_rank
# for sizex in 1000 2000 

for node in 1 
do
    for rank in 4 9 16 25 36 64
    do
        sizex=$(( 1000*$((node * rank)) ))
        echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi1d"
        srun -N $node --ntasks-per-node=$rank ./mpi2d -x $sizex -y 1000 -update 1 -t 128 -s 10
    done
done


for node in 2 
do
    for rank in 2 8 18 32 50
    do
        sizex=$(( 2000*$((node * rank)) ))
        echo "size = $sizex * 1000, freq = 1, steps = 128, seed = 10, node = $node, rank = $rank mpi1d"
        srun -N $node --ntasks-per-node=$rank ./mpi2d -x $sizex -y 1000 -update 1 -t 128 -s 10
    done
done