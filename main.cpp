#include <stdlib.h> 
#include <chrono>
#include <iostream>

#include "serial.cpp"


int main(int argc, char** argv) {

    int size = 10;
    int seed = 0;
    int steps = 100;
    int *data = new int[size * size];

    // random generate data
    for (int i = 0; i < size * size; ++i ) {
        data[i] = rand() % 2;
    }

    Board board(size);

    auto start_time = std::chrono::steady_clock::now();
    board.init_board(data, size);

    for (int timestamp = 0; timestamp < steps; ++timestamp ) {

    }

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    std::cout << "Simulation Time = " << seconds << " seconds." << std::endl;

    delete[] data;
}