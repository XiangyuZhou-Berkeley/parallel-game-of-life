#include <mpi.h>
#include <iostream>
#include <vector>
using namespace std;


int local_sizex = 0;
int local_sizey = 0;
int final_sizex = 0;
int update_frequency = 1;
int total_rank = 0;
//save data between row [begin_x,end_x]

int begin_x = 0;
int end_x;

vector<vector<int>> board;



//this data is sliced board data
void initiate(int rank, int sizex,int sizey, int* data, int ranks){
    local_sizey = sizey;
    local_sizex = sizex;
    total_rank = ranks;
    
    //calculate final_sizex
    //eg: local_size = 5, update = 2, final_x = 6
    //rank == 0, data saved [0, local_sizex -1]. bottom ghost area [local_sizex, final_sizex - 1]
    //rank == total_rank - 1, top ghost area [0,begin_x - 1].  data saved [beginx, final_sizex - 1].
    //else rank,top ghost area [0,begin_x - 1],data saved [beginx, beginx + local_sizex - 1].    .bottom ghost area [beginx + local_sizex, final_sizex - 1]
    if (rank == 0) {
        final_sizex = local_sizex + update_frequency;
        begin_x = 0;
        end_x = local_sizex - 1;

    } else if (rank == total_rank - 1){
        final_sizex = local_sizex + update_frequency;
        begin_x = update_frequency;
        end_x = final_sizex - 1;

    } else {
        begin_x = update_frequency;
        end_x = begin_x + local_sizex - 1;
        final_sizex = local_sizex + 2 * update_frequency;
    }
    
    //top_ghost area
    for (int i = 0; i < begin_x; ++i) {
        board.push_back(vector<int>(local_sizey, 0));
    }

    //put data into board
    for (int i = begin_x; i <= end_x; i++) {
        vector<int> temp_row(local_sizey, 0);
        for (int j = 0; j < local_sizey; j++) {
            temp_row[j] = data[(i - begin_x) * local_sizey + j];
        }
        board.push_back(temp_row);
    }

    //bottom_ghost area
    for (int i = end_x + 1; i < final_sizex; i++) {
        board.push_back(vector<int>(local_sizey, 0));
    }
}


void gather(int rank, int sizex, int sizey, int* data){
    int* recv_temp;
    int* send_data = new int[local_sizex * local_sizey];
    for (int i = begin_x; i <= end_x; i++) {
        for (int j = 0; j < sizey; j++) {
            send_data[(i - begin_x) * local_sizey + j] = board[i][j];
        }
    }
    MPI_Gather(send_data, local_sizex * local_sizey, MPI_INT, data, local_sizex * local_sizey, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] send_data;

}