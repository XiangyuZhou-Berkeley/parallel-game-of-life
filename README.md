# parallel-game-of-life

## Problem Description

According to [Wikipedia's article](https://en.wikipedia.org/wiki/Conway's_Game_of_Life): "The **Game of Life**, also known simply as **Life**, is a cellular automaton devised by the British mathematician John Horton Conway in 1970."

The board is made up of an `m x n` grid of cells, where each cell has an initial state: **live** (represented by a `1`) or **dead** (represented by a `0`). Each cell interacts with its [eight neighbors](https://en.wikipedia.org/wiki/Moore_neighborhood) (horizontal, vertical, diagonal) using the following four rules (taken from the above Wikipedia article):
<!-- self 1, neighbour <2 die
self 1, neighobour 2 or 3 , live
self 1 neighobour > 3 , die
self 0, neighbour 3  live -->
1. Any live cell with fewer than two live neighbors dies as if caused by under-population.
2. Any live cell with two or three live neighbors lives on to the next generation.
3. Any live cell with more than three live neighbors dies, as if by over-population.
4. Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction.

The next state is created by applying the above rules simultaneously to every cell in the current state, where births and deaths occur simultaneously. Given the current state of the `m x n` grid `board`, return *the next state*.



## Parallel Strategies:

### Using MPI to decompose grid into different processors



1. 1D decomposition
2. 2D decomposition

### Time independent update

instead of updating status at every time step, we are going to update every 2,3 or k steps.  So we could kind of get like time parallel, which could save the time of sending latency since the send amount of data is the same. But the computation time is differnt, so this method is a tradeoff between coputation and send latency. 

We will try to look up for the best k for both 1D decompostion and 2D decompostion.





# How to run our program:

 1.Create a build directory in mpi_1d, mpi_2d, serial

```
CMAKE ..
make

```

2. move according sh files from sh files directory to build directory.

3. ```
   sbatch -N 2 **.sh
   ```





# Map plot

Run ipynb in our plot directory

