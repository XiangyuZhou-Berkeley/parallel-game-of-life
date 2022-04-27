#include <vector>

using namespace std;

class Board{
public:
    Board(int size);
    void init_board(int *data, int size);
    
    private:
    int timestamp;
    const int size;
    vector<vector<int> > grid;
    void update();

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

int get_new_state(int i, int j) {
    
}

void Board::update() {
    vector<vector>> newgrid = grid;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++k) {

        } 
    }
}