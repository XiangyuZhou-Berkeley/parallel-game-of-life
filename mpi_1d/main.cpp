#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <random>
#include <cstring>

#include "common.h"


// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stoi(argv[iplace + 1]);
    }

    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}


int main(int argc, char** argv) {
    int steps = find_int_arg(argc, argv, "-t", 1000);
    int seed = find_int_arg(argc, argv, "-t", 10);
    int update_frequency = find_int_arg(argc, argv, "-update", 1);

    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int sizex = 10;
    int sizey = 10;


    int *data = new int[sizex * sizey];
    int *data_temp = new int[sizex * sizey];
    
    // random generate data as 0 or 1
    double p = 0.8; // the probability for generate as 1
    std::mt19937 gen(seed);
    srand(seed);
    std::discrete_distribution<> distrib({ 1-p, p });
                                        // ^^^  ^- probability for 1
                                        //  | probability for 0 
    // random generate data
    if (rank == 0){
        for (int i = 0; i < sizex * sizey; ++i ) {
            // data[i] = distrib(gen);
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

        //print input
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
    // int start_index = 0;
    // if (rank < residual) {
    //     start_index = rank * my_row * sizey;
    // } else{
    //     start_index = residual * (my_row+1) * sizey + (rank - residual + 1) * my_row * sizey;
    // }

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
    
    for (int timestamp = 0; timestamp < steps; timestamp += update_frequency ) {
        update(rank,timestamp);
    }

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
        for (int i = 0; i < sizex; ++i) {
            for (int j = 0; j < sizey; ++j){
               std::cout << data_temp[i * sizey + j] << " "; 
            }
            std::cout << std::endl;
        }
    }
    

    delete[] data;
    delete[] data_temp;
    MPI_Finalize();
}
