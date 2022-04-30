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
void initiate(int rank, int sizex,int sizey, int* data, int ranks, int frequency = 1){
    local_sizey = sizey;
    local_sizex = sizex;
    total_rank = ranks;
    update_frequency = frequency;
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


void update(int rank){
    /*send ghost area
    **recevie ghost area
    **decide which rows need to update and calculate 
    **calculate based on time
    */
    int num_transfer = update_frequency * local_sizey;
    int* s_upper_ghost = (int*) malloc((int)num_transfer  * sizeof(int));
    int* s_lower_ghost = (int*) malloc((int)num_transfer  * sizeof(int));
    int* r_upper_ghost = (int*) malloc((int)num_transfer  * sizeof(int));
    int* r_lower_ghost = (int*) malloc((int)num_transfer  * sizeof(int));



    MPI_Request requests[2];
    MPI_Status statuses[2];
    int num_send_row = min(update_frequency, local_sizex);
    int num_received_upper_ghost = 0;
    int num_received_lower_ghost = 0;
    int num_send_upper_ghost = num_send_row * local_sizey;
    int num_send_lower_ghost = num_send_row * local_sizey;
    //send ghost area
    
    //only s_lower_ghost
    if (rank == 0) {
        int num_temp = 0;
        for (int i = end_x - num_send_row + 1; i <= end_x; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_lower_ghost[num_temp++] = board[i][j];
            }
        }
    }

    if (rank == total_rank - 1) {
        //only s_upper_ghost,possible here only a few
        int num_temp = 0;
        for (int i = begin_x; i < begin_x + num_send_row; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_upper_ghost[num_temp++] = board[i][j];
            }
        }
    } 

    if (rank > 0 && rank < total_rank - 1) {
        //lower ghost
        int num_temp = 0;
        for (int i = end_x - num_send_row + 1; i <= end_x; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_lower_ghost[num_temp++] = board[i][j];
            }
        }

        //upper ghost
        num_temp = 0;
        for (int i = begin_x; i < begin_x + num_send_row; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_upper_ghost[num_temp++] = board[i][j];
            }
        }

        //put into the grid
        
    }
    

    if (rank > 0) {
        MPI_Isend(s_upper_ghost, num_send_upper_ghost, MPI_INT, rank-1, rank, MPI_COMM_WORLD, &requests[0]);
        // receive upper ghost from previous
        MPI_Probe(rank-1 , rank - 1 , MPI_COMM_WORLD , &statuses[0]);
        MPI_Get_count(&statuses[0], MPI_INT, &num_received_upper_ghost);
        MPI_Recv(r_upper_ghost,num_received_upper_ghost, MPI_INT, rank-1, rank-1, MPI_COMM_WORLD, &statuses[0]);
    }

    if (rank < total_rank - 1){
        MPI_Isend(s_lower_ghost, num_send_lower_ghost, MPI_INT, rank+1, rank, MPI_COMM_WORLD, &requests[1]);
        // receive lower ghost from followgin processor
        MPI_Probe(rank+1 , rank+1 , MPI_COMM_WORLD , &statuses[1]);
        MPI_Get_count(&statuses[1], MPI_INT, &num_received_lower_ghost);
        MPI_Recv(r_lower_ghost, num_received_lower_ghost, MPI_INT, rank+1, rank+1, MPI_COMM_WORLD, &statuses[1]);
    }
}


void calculate_all(int rank) {
    int begin_calculate_x = 0;
    int end_calculate_x = 0;
    if (rank > 0) {
        begin_calculate_x = 1;
    } else {
        begin_calculate_x = begin_x;
    }
    if (rank < total_rank - 1) {
        end_calculate_x = final_sizex - 1 - 1;
    } else {
        end_calculate_x = end_x;
    }

    for (time = 0; time < update_frequency; time++) {
        
    }
}

void calculate_once(){
    
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