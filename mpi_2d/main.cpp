#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <fstream>
#include <string>

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

void output_to_file(std::ofstream& out, int *data, int sizex, int sizey) {
    for (int i = 0; i < sizex; ++i) {
        for (int j = 0; j < sizey; ++j) {
            out << data[i * sizey + j] << " ";
        }
        out << std::endl;
    }
}



int main(int argc, char** argv) {
    int steps = find_int_arg(argc, argv, "-t", 1000);
    int update_frequency = find_int_arg(argc, argv, "-update", 1);
    int seed = find_int_arg(argc, argv, "-s", 10);
    char* filename = find_string_option(argc, argv, "-i", nullptr);
    char* outputfile = find_string_option(argc, argv, "-o", nullptr);
    
    int sizex = 10;
    int sizey = 10;

    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *data;
    int *data_temp;
    int *final_output;

    if (filename == nullptr) {
        // cout << "since no file found, we use random generate" << endl;
        // random generate data as 0 or 1
        srand(seed);
        sizex = find_int_arg(argc, argv, "-x", 10);
        sizey = find_int_arg(argc, argv, "-y", 10);
        data = new int[sizex * sizey];
        data_temp = new int[sizex * sizey];
        final_output = new int[sizex * sizey];

        double p = 0.8; // the probability for generate as 1
        std::mt19937 gen(seed);
        std::discrete_distribution<> distrib({ 1-p, p });
                                            // ^^^  ^- probability for 1
                                            //  | probability for 0 
        if (rank == 0){
            for (int i = 0; i < sizex * sizey; ++i ) {
                // data[i] = distrib(gen);
                data[i] = rand() % 2;
            }

            //  print input
            // for (int i = 0; i < sizex; ++i) {
            //     for (int j = 0; j < sizey; ++j){
            //         std::cout << data[i * sizey + j] << " "; 
            //     }
            //     std::cout << std::endl;
            // }
        }
    } else {
        // generate from file initialize
        // we need to change sizex, sizey according to the size of the file
        sizex = 0;
        sizey = 0; 
        ifstream readfile_size(filename);
        if ( readfile_size.is_open() ){
            string fileline_size;  
            while (getline(readfile_size,fileline_size))
            {   
                if (fileline_size.length() > sizey){
                    sizey = fileline_size.length();
                }
                sizex = sizex + 1;
            }
        } else {
            cout << "No such file to read size, try again." << endl;
        }
        sizex = (sizex / 100 + 1) * 100;
        sizey = (sizey / 100 + 1) * 100;
        data = new int[sizex * sizey];
        data_temp = new int[sizex * sizey];
        final_output = new int[sizex * sizey];
        if (rank == 0){
            // first initialize to make usre every grid has a value
            for (int i = 0; i < sizex * sizey; ++i ) {
                data[i] = 0;
            }

            //string filename = "../../initialize_data/bhept1.txt";
            ifstream readfile(filename);
            if ( readfile.is_open() ){
                string fileline;
                int row = 0;
                int column = 0;  
                char item;
                while (getline(readfile,fileline))
                {
                    for (int i = 0; i < fileline.length(); ++i){
                        item = fileline[i];
                        column = column + 1;
                        if (item == '*'){
                            data[row * sizey + column] = 1;
                        } 
                    }
                    row = row + 1;
                    column = 0;
                }

                //print input
                // for (int i = 0; i < sizex; ++i) {
                //     for (int j = 0; j < sizey; ++j){
                //         std::cout << data[i * sizey + j] << " "; 
                //     }
                //     std::cout << std::endl;
                // }
            } else {
                cout << "No such file, try again." << endl;
            }
        }
    }


    //broadcast data to every processor
    MPI_Bcast(data, sizex * sizey, MPI_INT, 0, MPI_COMM_WORLD);
    int proc_per_row = sqrt(num_procs);
    int proc_per_col = sqrt(num_procs);
    // if (rank == 0) {
    //     std::cout << proc_per_row << std::endl;
    // }
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

    auto start_time = std::chrono::steady_clock::now();
   

    if (rank < total_rank){
        initiate(rank, my_row, my_col,reshaped_data[rank].data(), total_rank,update_frequency);
        for (int i = 0; i < steps; i += update_frequency) {
            update(rank, i, proc_per_row);
        }
        gather(rank,data_temp,displacement,recvcounts);
    }

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
        auto end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;
        double seconds = diff.count();
        std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

        // // print output
        // if (outputfile != nullptr) {
        //     std::ofstream out(outputfile);
        //     output_to_file(out, final_output, sizex, sizey);
        // }
    }
    

    delete[] data;
    delete[] data_temp;
    delete[] final_output;
    MPI_Finalize();
}