#include <stdlib.h> 
#include <chrono>
#include <iostream>

#include "serial.cpp"


int main(int argc, char** argv) {

    
    int sizex = 5;
    int sizey = 5;
    int seed = 0;
    int steps = 2;
    int *data = new int[sizex * sizey];

    // random generate data
    // for (int i = 0; i < size * size; ++i ) {
    //     data[i] = rand() % 2;
    // }
    int *temp = new int[sizex * sizey] {0,1,0, 0, 1,0,0,1, 1, 0,1,1,1, 0 ,1, 1,0, 0, 0, 1,0,1,0,0,1};

    data = temp;

    Board board(sizex, sizey);

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, sizex, sizey);

    for (int timestamp = 0; timestamp < steps; ++timestamp ) {
        board.update();
    }
    std::cout << "Finished simulation" << std::endl;
    board.print_board();

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

    delete[] data;
}