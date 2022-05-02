#include <stdlib.h> 
#include <chrono>
#include <iostream>
#include <random>

#include "serial.cpp"


int main(int argc, char** argv) {

    
    int sizex = 10;
    int sizey = 10;
    int seed = 10;
    int steps = 100;
    int *data = new int[sizex * sizey];

    // random generate data as 0 or 1
    double p = 0.8; // the probability for generate as 1
    std::mt19937 gen(seed);
    std::discrete_distribution<> distrib({ 1-p, p });
                                        // ^^^  ^- probability for 1
                                        //  | probability for 0 
    for (int i = 0; i < sizex * sizey; ++i ) {
        data[i] = distrib(gen);
    }
    //int *temp = new int[sizex * sizey] {0,1,0, 0, 1,0,0,1, 1, 0,1,1,1, 0 ,1, 1,0, 0, 0, 1,0,1,0,0,1};

    //data = temp;
    
    Board board(sizex, sizey);
    

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, sizex, sizey);
    board.print_board();
    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        board.update();
    }
    
    board.print_board();

    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Finished simulation" << std::endl;
    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

    delete[] data;
}