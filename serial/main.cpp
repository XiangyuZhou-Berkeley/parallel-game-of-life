#include <stdlib.h> 
#include <chrono>
#include <iostream>

#include "serial.cpp"


int main(int argc, char** argv) {

    
    int sizex = 1000;
    int sizey = 1000;
    int seed = 10;
    int steps = 100;
    int *data = new int[sizex * sizey];
    srand(seed);
    // random generate data
    for (int i = 0; i < sizex * sizey; ++i ) {
        data[i] = rand() % 2;
    }
    //int *temp = new int[sizex * sizey] {0,1,0, 0, 1,0,0,1, 1, 0,1,1,1, 0 ,1, 1,0, 0, 0, 1,0,1,0,0,1};

    //data = temp;
    
    Board board(sizex, sizey);
    

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, sizex, sizey);
    //board.print_board();
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