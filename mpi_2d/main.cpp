#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <cstring>

#include "common.h"
using namespace std;

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
    int update_frequency = find_int_arg(argc, argv, "-update", 1);
    int seed = find_int_arg(argc, argv, "-s", 10);

    int sizex = find_int_arg(argc, argv, "-x", 10);
    int sizey = find_int_arg(argc, argv, "-y", 10);


    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *data = new int[sizex * sizey];
    int *data_temp = new int[sizex * sizey];
    int *final_output = new int[sizex * sizey];
    srand(seed);
    // random generate data
    if (rank == 0){
        for (int i = 0; i < sizex * sizey; ++i ) {
            data[i] = rand() % 2;
        }

        //  print input
        for (int i = 0; i < sizex; ++i) {
            for (int j = 0; j < sizey; ++j){
                std::cout << data[i * sizey + j] << " "; 
            }
            std::cout << std::endl;
        }
    }


    //broadcast data to every processor
    MPI_Bcast(data, sizex * sizey, MPI_INT, 0, MPI_COMM_WORLD);
    int proc_per_row = sqrt(num_procs);
    int proc_per_col = sqrt(num_procs);
    int total_rank = proc_per_col * proc_per_col;
    int row_per_proc = sizex / proc_per_row;
    int col_per_proc = sizey / proc_per_row;
    int row_residual = sizex - row_per_proc * proc_per_row;
    int col_residual = sizey - col_per_proc * proc_per_row; 
    //here think about the edge case like we have 5 rows, two processors, row_per_proc = 3, but processor 1 only have 2 rows
    
    int my_row = row_per_proc;
    int my_col = col_per_proc;
    int rank_row = rank / proc_per_col;
    int rank_col = rank % proc_per_col;

    //we will only use sqrt(num_procs) * sqrt(num_procs) processors
    if (rank_row < row_residual){
        my_row = my_row + 1;
    }
    if (rank_col < col_residual) {
        my_col = my_col + 1;
    }
    //we only use these ranks,get which line and column to begin for initiate and gather 
    int* row_displacement = new int[row_per_proc];
    int* col_displacement = new int[col_per_proc];
    for (int i = 0; i < row_per_proc; i++) {
        if (i < row_residual) {
            row_displacement[i] = i * (row_per_proc + 1);
        } else {
            row_displacement[i] = row_residual * (row_per_proc + 1) + (i - row_residual) * row_per_proc;
        }
        // if (rank == 0) {
        // }
    }

    for (int i = 0; i < col_per_proc; i++) {
        if (i < col_residual) {
            col_displacement[i] = i * (col_per_proc + 1);
        } else {
            col_displacement[i] = col_residual * (row_per_proc + 1) + (i - col_residual) * row_per_proc;
        }
        // if (rank == 0) {
        //     std::cout << col_displacement[i] <<std::endl;
        // }
    }

    //every vector is savxes 1d data for every processor
    vector<vector<int>> reshaped_data;
    //recvcount and siplacement for gather
    int * recvcounts = new int[total_rank];
    int* displacement = new int[total_rank];

    for (int i = 0; i < total_rank; i++) {
        int rank_row_temp = i / proc_per_col;
        int rank_col_temp = i % proc_per_col;
        
        //[start_row,end_row)
        int start_row = row_displacement[rank_row_temp];
        int end_row;
        int start_col = col_displacement[rank_col_temp];
        int end_col;
        if (rank_row_temp == proc_per_row - 1) {
            end_row = sizex;
        } else {
            end_row = row_displacement[rank_row_temp + 1];
        }

        if (rank_col_temp == proc_per_col - 1) {
            end_col = sizey;
        } else {
            end_col = col_displacement[rank_col_temp + 1];
        }
        int size_temp = (end_row - start_row) * (end_col - start_col);
        recvcounts[i] = size_temp;
        // if (rank == 0){
        //     std::cout<< i << " " << size_temp << std::endl;
        // }
        
        vector<int> temp_rank_data;
        for (int r = start_row; r < end_row; r++) {

            for (int c = start_col; c < end_col; c++) {
                temp_rank_data.push_back(data[r * sizey + c]);
            }
        }
        reshaped_data.push_back(temp_rank_data);
    }
    displacement[0] = 0;
    for (int i = 1; i < total_rank; i++) {
        displacement[i] = displacement[i - 1] + recvcounts[i - 1];
        // if (rank == 0) {
        //     cout<< i <<" number: " << recvcounts[i] << " displacement" << displacement[i] << endl;
        // }
    }
    
    // debug to test reshaped data
    // if (rank == 0) {
    //     for (int i = 0; i < reshaped_data.size(); i++) {
    //         std::cout << "rank" << i << ":" << std::endl;
    //         for (int j = 0; j < reshaped_data[i].size(); j++) {
    //             std::cout << reshaped_data[i][j] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // if (rank < proc_per_row * proc_per_col) {
        
    // }




    
    
    

    // int* displacement = new int[num_procs];
    // int * recvcounts = new int[num_procs];
    // for (int i = 0; i < num_procs; i++) {
    //     if (i < residual) {
    //         displacement[i] = i * (row_per_proc + 1) * sizey;
    //         recvcounts[i] = (row_per_proc + 1) * sizey;
    //     } else {
    //         displacement[i] = (residual * (row_per_proc + 1) + (i - residual) * row_per_proc) * sizey;
    //         recvcounts[i] = row_per_proc * sizey;
    //     }
    //     if (rank == 0) {
    //     }
    // }

    if (rank < total_rank){
        initiate(rank, my_row, my_col,reshaped_data[rank].data(), total_rank,update_frequency);
        for (int i = 0; i < steps; i += update_frequency) {
            update(rank, i, proc_per_row);
        }
        gather(rank,data_temp,displacement,recvcounts);
        auto start_time = std::chrono::steady_clock::now();
    }
    
    //debug to make sure gather is right:
    // if (rank == 0) {
    //     std::cout << endl;
    //     for (int i = 0; i < total_rank; i++) {
    //         for (int j = 0; j < recvcounts[i]; j++) {
    //             std::cout << data_temp[displacement[i] + j] <<" ";
    //         }
    //         std::cout << endl;
    //     }
    // }

    // for (int timestamp = 0; timestamp < steps; ++timestamp ) {
    //     update(rank,timestamp);
    // }
    
    // gather(rank, sizex, sizey, data_temp, displacement, recvcounts);

    // auto end_time = std::chrono::steady_clock::now();
    // if (rank == 0){
    //     std::cout << "Finished simulation" << std::endl;
    // }


    //reput data into 2d grid like
    if(rank == 0){
        for (int i = 0; i < total_rank; i++) {
            int rank_row_temp = i / proc_per_col;
            int rank_col_temp = i % proc_per_col;
            
            //[start_row,end_row)
            int start_row = row_displacement[rank_row_temp];
            int end_row;
            int start_col = col_displacement[rank_col_temp];
            int end_col;
            if (rank_row_temp == proc_per_row - 1) {
                end_row = sizex;
            } else {
                end_row = row_displacement[rank_row_temp + 1];
            }

            if (rank_col_temp == proc_per_col - 1) {
                end_col = sizey;
            } else {
                end_col = col_displacement[rank_col_temp + 1];
            }
            int temp_index = 0;
            for (int r = start_row; r < end_row; r++) {
                for (int c = start_col; c < end_col; c++) {
                    final_output[r * sizey + c] = data_temp[displacement[i] + temp_index];
                    temp_index++;
                }
            }
        }
    }



    if (rank == 0){
        // std::chrono::duration<double> diff = end_time - start_time;
        // double seconds = diff.count();
        // std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

        // print output
        std::cout << std::endl;
        for (int i = 0; i < sizex; ++i) {
            for (int j = 0; j < sizey; ++j){
               std::cout << final_output[i * sizey + j] << " "; 
            }
            std::cout << std::endl;
        }
    }
    

    delete[] data;
    delete[] data_temp;
    delete[] final_output;
    MPI_Finalize();
}