#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include <string>

#include "serial.cpp"


int main(int argc, char** argv) {

    
    int sizex = 10;
    int sizey = 10;
    int seed = 10;
    int steps = 100;

    // generate from file initialize
    sizex = 300; // this is fixed
    sizey = 300; // this is fixed
    // // the following size is just for test, use test_jxy.txt instead:
    // sizex = 10;
    // sizey = 10;

    int *data = new int[sizex * sizey];

    // // random generate data as 0 or 1
    // double p = 0.8; // the probability for generate as 1
    // std::mt19937 gen(seed);
    // std::discrete_distribution<> distrib({ 1-p, p });
    //                                     // ^^^  ^- probability for 1
    //                                     //  | probability for 0 
    // for (int i = 0; i < sizex * sizey; ++i ) {
    //     data[i] = distrib(gen);
    // }
    
    
    // generate from file
    // first initialize to make usre every grid has a value
    for (int i = 0; i < sizex * sizey; ++i ) {
        data[i] = 0;
    }

    string filename = "../../initialize_data/bhept1.txt";
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

    
    Board board(sizex, sizey);
    

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, sizex, sizey);
    board.print_board();
    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        board.update();
    }
    
    //board.print_board();

    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Finished simulation" << std::endl;
    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

    delete[] data;
}