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
    int sizex = 1000;
    int sizey = 1000;
    int seed = 10;
    int steps = 100;
    int update_frequency = 1;
    int *data = new int[sizex * sizey];
    int *data_temp = new int[sizex * sizey];
    srand(seed);
    // random generate data
    if (rank == 0){
        for (int i = 0; i < sizex * sizey; ++i ) {
            data[i] = rand() % 2;
        }
    }


    //broadcast data to every processor
    MPI_Bcast(data, sizex * sizey, MPI_INT, 0, MPI_COMM_WORLD);
    int proc_per_row = sqrt(num_procs);
    int proc_per_col = sqrt(num_procs);
    int row_per_proc = sizex / proc_per_row;
    int col_per_proc = sizey / proc_per_row;
    int row_residual = sizex - row_per_proc * num_procs;
    int col_residual = sizey - col_per_proc * num_procs; 
    //here think about the edge case like we have 5 rows, two processors, row_per_proc = 3, but processor 1 only have 2 rows
    
    int my_row = row_per_proc;
    int my_col = col_per_proc;
    int rank_row = rank / proc_per_col;
    int rank_col = rank % proc_per_col;

    //we will only use sqrt(num_procs) * sqrt(num_procs) processors
    if (rank_row < row_residual)
    if (rank < residual) {
        my_row = my_row + 1;
    }
    

    int start_index = 0;
    if (rank < residual) {
        start_index = rank * (row_per_proc + 1) * sizey;
    } else {
        start_index = (residual * (row_per_proc + 1) + (rank - residual) * row_per_proc) * sizey;
    }
    
    

    int* displacement = new int[num_procs];
    int * recvcounts = new int[num_procs];
    for (int i = 0; i < num_procs; i++) {
        if (i < residual) {
            displacement[i] = i * (row_per_proc + 1) * sizey;
            recvcounts[i] = (row_per_proc + 1) * sizey;
        } else {
            displacement[i] = (residual * (row_per_proc + 1) + (i - residual) * row_per_proc) * sizey;
            recvcounts[i] = row_per_proc * sizey;
        }
        if (rank == 0) {
            //std::cout << displacement[i] <<std::endl;
           // std::cout << recvcounts[i] <<std::endl;
        }
    }
    auto start_time = std::chrono::steady_clock::now();
    initiate(rank, my_row, sizey, data + start_index, num_procs,update_frequency);
    


    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        update(rank,timestamp);
    }
    
    // board.print_board();
    //update(rank);

    // std::cout << "Rank " << rank << " finished update" << std::endl;
    gather(rank, sizex, sizey, data_temp, displacement, recvcounts);

    auto end_time = std::chrono::steady_clock::now();
    if (rank == 0){
        std::cout << "Finished simulation" << std::endl;
    }
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

        // print output
        // for (int i = 0; i < sizex; ++i) {
        //     for (int j = 0; j < sizey; ++j){
        //        std::cout << data_temp[i * sizey + j] << " "; 
        //     }
        //     std::cout << std::endl;
        // }
    }
    

    delete[] data;
    MPI_Finalize();
}