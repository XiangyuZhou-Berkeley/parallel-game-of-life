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
    int sizex = 5;
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
    int row_per_proc = floor((sizex - 1) / num_procs) + 1;
    //here think about the edge case like we have 5 rows, two processors, row_per_proc = 3, but processor 1 only have 2 rows
    
    int my_row = row_per_proc;
    if (rank == num_procs - 1) {
        my_row = sizex - (num_procs - 1) * row_per_proc;
    }
    initiate(rank, my_row, sizey, data + rank * row_per_proc * sizey , num_procs,1);
    auto start_time = std::chrono::steady_clock::now();


    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        // board.update();
    }
    if (rank == 0){
        std::cout << "Finished simulation" << std::endl;
    }
    
    // board.print_board();
    gather(rank, sizex, sizey, data_temp);
    auto end_time = std::chrono::steady_clock::now();
    
    if (rank == 0) {
         for (int i = 0; i < sizex; ++i) {
            for (int j = 0; j < sizey; ++j){
                std::cout << data_temp[i * sizey + j] << " "; 
            }
            std::cout << std::endl;
        }
    }


    if (rank == 0){
        std::chrono::duration<double> diff = end_time - start_time;
        double seconds = diff.count();
        std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;
    }
    

    delete[] data;
    MPI_Finalize();
}