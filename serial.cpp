#include <vector>

using namespace std;

class Board{
public:
    Board(int size);
    void init_board(int *data, int size);
    
    private:
    int timestamp;
    const int size;
    vector<vector<int>> grid;
    void update();
    int get_new_state(int i, int j);

};


Board::Board(int size): timestamp(0), size(size), grid(size, vector<int>(size)) {

}

void Board::init_board(int *data, int size){
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            grid[i][j] = data[size*i+j];
        }
    }
}

int Board::get_new_state(int i, int j) {
    return grid[i][j];
}

void Board::update() {

    vector<vector<int>> newgrid = grid;

    for (int i = 0; i < size; ++i) {
        int row_start = i - 1 >= 0 ? i - 1 : 0;
        int row_end = i + 1 < size ? i + 1 : size - 1;
        for (int j = 0; j < size; ++j) {
            int col_start = j - 1 >= 0 ? j - 1 : 0;
            int col_end = j + 1 < size ? j + 1 : size - 1;
            int alive_neighbour = 0;
            //inlcude itself
            for (int row = row_start; row < row_end; row++) {
                for (int col = col_start; col < col_end; col++) {
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