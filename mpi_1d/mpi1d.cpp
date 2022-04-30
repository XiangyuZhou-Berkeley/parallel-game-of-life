#include <mpi.h>
#include <iostream>
#include <vector>
using namespace std;


int local_sizex = 0;
int local_sizey = 0;
int final_sizex = 0;
int update_frequency = 1;
int total_rank = 0;


int begin_x = 0;
int end_x;
vector<vector<int>> upper_ghost;
vector<vector<int>> lower_ghost;

//only save[0,local_sizex)
vector<vector<int>> board;



void initiate(int rank, int sizex,int sizey, int* data, int ranks, int frequency = 1){
    local_sizey = sizey;
    local_sizex = sizex;
    total_rank = ranks;
    update_frequency = frequency;
   

    //put data into board
    for (int i = 0; i < local_sizex; i++) {
        vector<int> temp_row(local_sizey, 0);
        for (int j = 0; j < local_sizey; j++) {
            temp_row[j] = data[i * local_sizey + j];
        }
        board.push_back(temp_row);
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
        for (int i = local_sizex - num_send_row; i < local_sizex; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_lower_ghost[num_temp++] = board[i][j];
            }
        }
    }

    if (rank == total_rank - 1) {
        //only s_upper_ghost,possible here only a few
        int num_temp = 0;
        for (int i = 0; i < num_send_row; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_upper_ghost[num_temp++] = board[i][j];
            }
        }
    } 

    if (rank > 0 && rank < total_rank - 1) {
        //lower ghost
        int num_temp = 0;
        for (int i = local_sizex - num_send_row; i < local_sizex; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_lower_ghost[num_temp++] = board[i][j];
            }
        }

        //upper ghost
        num_temp = 0;
        for (int i = 0; i < num_send_row; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_upper_ghost[num_temp++] = board[i][j];
            }
        }


    }
    
    //send recv and put into grid
    if (rank > 0) {
        MPI_Isend(s_upper_ghost, num_send_upper_ghost, MPI_INT, rank-1, rank, MPI_COMM_WORLD, &requests[0]);
        // receive upper ghost from previous
        MPI_Probe(rank-1 , rank - 1 , MPI_COMM_WORLD , &statuses[0]);
        MPI_Get_count(&statuses[0], MPI_INT, &num_received_upper_ghost);
        MPI_Recv(r_upper_ghost,num_received_upper_ghost, MPI_INT, rank-1, rank-1, MPI_COMM_WORLD, &statuses[0]);
        //put in first receive_row
        int receive_row = num_received_upper_ghost / local_sizey;
        int num_temp = 0;
        for (int i = 0; i < receive_row; i++) {
            vector<int> temp_row(local_sizey, 0);
            for (int j = 0; j < local_sizey; j++) {
                temp_row[j] = r_upper_ghost[i * local_sizey + j];
            }
            upper_ghost.push_back(temp_row);
        }
    }

    if (rank < total_rank - 1){
        MPI_Isend(s_lower_ghost, num_send_lower_ghost, MPI_INT, rank+1, rank, MPI_COMM_WORLD, &requests[1]);
        // receive lower ghost from followgin processor
        MPI_Probe(rank+1 , rank+1 , MPI_COMM_WORLD , &statuses[1]);
        MPI_Get_count(&statuses[1], MPI_INT, &num_received_lower_ghost);
        MPI_Recv(r_lower_ghost, num_received_lower_ghost, MPI_INT, rank+1, rank+1, MPI_COMM_WORLD, &statuses[1]);
        int receive_row = num_received_lower_ghost / local_sizey;
        int num_temp = 0;
        //put in last receive_row 
        for (int i = 0; i < receive_row; i++) {
            vector<int> temp_row(local_sizey, 0);
            for (int j = 0; j < local_sizey; j++) {
                temp_row[j] = r_lower_ghost[i * local_sizey + j];
            }
            lower_ghost.push_back(temp_row);
        }
    }
}


void calculate_all(int rank) {
    
}



void gather(int rank, int sizex, int sizey, int* data){
    int* recv_temp;
    int* send_data = new int[local_sizex * local_sizey];
    for (int i = 0; i < local_sizex; i++) {
        for (int j = 0; j < sizey; j++) {
            send_data[i * local_sizey + j] = board[i][j];
        }
    }
    MPI_Gather(send_data, local_sizex * local_sizey, MPI_INT, data, local_sizex * local_sizey, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] send_data;

}