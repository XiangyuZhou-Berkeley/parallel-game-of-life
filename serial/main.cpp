#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <cstring>

#include "serial.cpp"

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
    int seed = find_int_arg(argc, argv, "-s", 10);
    int sizex = 10;
    int sizey = 10;

    char* filename = find_string_option(argc, argv, "-i", nullptr);
    char* outputfile = find_string_option(argc, argv, "-o", nullptr);

    int *data;

    if (filename == nullptr) {
        cout << "since no file found, we use random generate" << endl;
        srand(seed);

        sizex = find_int_arg(argc, argv, "-x", 10);
        sizey = find_int_arg(argc, argv, "-y", 10);

        data = new int[sizex * sizey];

        // random generate data as 0 or 1
        double p = 0.8; // the probability for generate as 1
        std::mt19937 gen(seed);
        std::discrete_distribution<> distrib({ 1-p, p });
                                            // ^^^  ^- probability for 1
                                            //  | probability for 0 
        for (int i = 0; i < sizex * sizey; ++i ) {
            // data[i] = distrib(gen);
            data[i] = rand() % 2;
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
        } else {
            cout << "No such file, try again." << endl;
        }
    }

    
    Board board(sizex, sizey);
    

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, sizex, sizey);
    board.print_board();
    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        board.update();
    }
    
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Finished simulation" << std::endl;

    if (outputfile != nullptr) {
        std::ofstream out(outputfile);
        board.output_to_file(out);
    }
    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

    delete[] data;
}