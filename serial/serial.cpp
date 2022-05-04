#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Board{
public:
    Board(int sizex, int sizey);
    void init_board(int *data, int sizex, int sizey);
    void update();
    void print_board();
    void output_to_file(std::ofstream& out);
private:
    int timestamp;
    const int sizex; // row
    const int sizey; // col
    vector<vector<int>> grid;
    int get_new_state(int i, int j);

};


Board::Board(int sizex, int sizey): timestamp(0), sizex(sizex), sizey(sizey),
                        grid(sizex, vector<int>(sizey)) {

}

void Board::init_board(int *data, int sizex, int sizey){
    for (int i = 0; i < sizex; ++i) {
        for (int j = 0; j < sizey; ++j) {
            grid[i][j] = data[sizey*i+j];
        }
    }
}

int Board::get_new_state(int i, int j) {
    return grid[i][j];
}

void Board::update() {

    vector<vector<int>> newgrid = grid;

    for (int i = 0; i < sizex; ++i) {
        int row_start = i - 1 >= 0 ? i - 1 : 0;
        int row_end = i + 1 < sizex ? i + 1 : sizex - 1;
        for (int j = 0; j < sizey; ++j) {
            int col_start = j - 1 >= 0 ? j - 1 : 0;
            int col_end = j + 1 < sizey ? j + 1 : sizey - 1;
            int alive_neighbour = 0;
            //inlcude itself
            for (int row = row_start; row <= row_end; row++) {
                for (int col = col_start; col <= col_end; col++) {
                    alive_neighbour = alive_neighbour + grid[row][col];
                }
            }
            alive_neighbour = alive_neighbour - grid[i][j];
            
            //logistic 
            //self is 1
            if (grid[i][j] == 1) {
                if (alive_neighbour == 2 || alive_neighbour == 3){
                } else {
                    newgrid[i][j] = 0;
                }
            } else {
                if (alive_neighbour == 3) {
                    newgrid[i][j] = 1;
                }
            }
        } 
    }
    grid = newgrid;
}

void Board::print_board(){
    for (int i = 0; i < sizex; ++i) {
        for (int j = 0 ; j < sizey; ++j) {
            cout << grid[i][j] << " ";
        }
        cout << endl;
    }

}

void Board::output_to_file(std::ofstream& out) {
    for (int i = 0; i < sizex; ++i) {
        for (int j = 0; j < sizey; ++j) {
            out << grid[i][j] << " ";
        }
        out << std::endl;
    }
}