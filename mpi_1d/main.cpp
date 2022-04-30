#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include "common.h"
// #include "mpi1d.cpp"


int main(int argc, char** argv) {
    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int sizex = 7;
    int sizey = 5;
    int seed = 0;
    int steps = 2;
    int *data = new int[sizex * sizey];
    int *data_temp = new int[sizex * sizey];
    // random generate data
    if (rank == 0){
        for (int i = 0; i < sizex * sizey; ++i ) {
            data[i] = rand() % 2;
        }

        //manual test case
        // 1 0 0 1 0 
        // 0 0 0 0 0 
        // 1 1 1 1 0 
        // 0 0 0 0 0 
        // 0 0 0 0 1
        // int *temp = new int [sizex*sizey]{1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1};
        // data = temp;
        for (int i = 0; i < sizex; ++i) {
            for (int j = 0; j < sizey; ++j){
                std::cout << data[i * sizey + j] << " "; 
            }
            std::cout << std::endl;
        }
    }

    // int *temp = new int[sizex * sizey] {0,1,0, 0, 1,0,0,1, 1, 0,1,1,1, 0 ,1, 1,0, 0, 0, 1,0,1,0,0,1};

    // data = temp;

    //broadcast data to every processor
    MPI_Bcast(data, sizex * sizey, MPI_INT, 0, MPI_COMM_WORLD);
    int row_per_proc = sizex / num_procs;
    int residual = sizex - row_per_proc * num_procs;
    //here think about the edge case like we have 5 rows, two processors, row_per_proc = 3, but processor 1 only have 2 rows
    
    int my_row = row_per_proc;
    if (rank < residual) {
        my_row = my_row + 1;
    }
    int start_index = 0;
    if (rank < residual) {
        start_index = rank * my_row * sizey;
    } else{
        start_index = residual * my_row * sizey + (rank - residual + 1) * my_row * sizey;
    }
    int* displacement = new int[num_procs];
    for (int i = 0; i < num_procs; i++) {
        if (i < residual) {
            displacement[i] = i * (row_per_proc + 1);
        } else {
            displacement[i] = residual * (row_per_proc + 1) + (i - residual) * row_per_proc;
        }
        if (rank == 0) {
            std::cout << displacement[i] <<std::endl;
        }
    }
    
    initiate(rank, my_row, sizey, data + start_index, num_procs,1);
    auto start_time = std::chrono::steady_clock::now();


    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        // board.update();
    }
    if (rank == 0){
        std::cout << "Finished simulation" << std::endl;
    }
    
    // board.print_board();
    update(rank);

    std::cout << "Rank " << rank << " finished update" << std::endl;
    gather(rank, sizex, sizey, data_temp);

    auto end_time = std::chrono::steady_clock::now();
    
    // if (rank == 0) {
    //      for (int i = 0; i < sizex; ++i) {
    //         for (int j = 0; j < sizey; ++j){
    //             std::cout << data_temp[i * sizey + j] << " "; 
    //         }
    //         std::cout << std::endl;
    //     }
    // }


    if (rank == 0){
        std::chrono::duration<double> diff = end_time - start_time;
        double seconds = diff.count();
        std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;
    }
    

    delete[] data;
    MPI_Finalize();
}