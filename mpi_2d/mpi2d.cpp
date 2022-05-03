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


//only save[0,local_sizex)
vector<vector<int>> board;



void initiate(int rank, int sizex,int sizey, int* data, int ranks, int frequency = 1){
    local_sizey = sizey;
    local_sizex = sizex;
    total_rank = ranks;
    update_frequency = frequency;

    // std::cout << "[Initiate]： rank " << rank << "local_sizex: " << local_sizex << std::endl;
   

    //put data into board
    for (int i = 0; i < local_sizex; i++) {
        vector<int> temp_row(local_sizey, 0);
        for (int j = 0; j < local_sizey; j++) {
            temp_row[j] = data[i * local_sizey + j];
        }
        board.push_back(temp_row);
    }

}


void calculate_all(int rank, int step, vector<vector<int>> &upper_ghost, vector<vector<int>> &lower_ghost,vector<vector<int>> left_ghost,
    vector<vector<int>> right_ghost,
    vector<vector<int>> upper_left_ghost,
    vector<vector<int>> upper_right_ghost,
    vector<vector<int>> lower_left_ghost,
    vector<vector<int>> lower_right_ghost) {
    int upper_x = upper_ghost.size();
    int lower_x = lower_ghost.size();
    int new_x = upper_x + local_sizex + lower_x;
    int left_y = left_ghost.size() == 0 ? 0 : left_ghost[0].size();
    int right_y = right_ghost.size() == 0 ? 0 : right_ghost[0].size();
    int new_y = left_y + local_sizey + right_y;

    vector<vector<int> > new_board (new_x, vector<int>(new_y, 0));

    for (int i = 0; i < upper_left_ghost.size(); ++i) {
        for (int j = 0; j < upper_left_ghost[i].size(); ++j) {
            new_board[i][j] = upper_left_ghost[i][j];
        }
    }

    for (int i = 0; i < upper_ghost.size(); ++i) {
        for (int j = 0; j < upper_ghost[i].size(); ++j) {
            new_board[i][j+left_y] = upper_ghost[i][j];
        }
    }

    // upper right
    for (int i = 0; i < upper_right_ghost.size(); ++i ){
        for (int j = 0; j < upper_right_ghost[i].size(); ++j) {
            new_board[i][j+left_y+local_sizey] = upper_left_ghost[i][j];
        }
    }

    // left 
    for (int i = 0; i < left_ghost.size(); ++i) {
        for (int j = 0; j < left_ghost[i].size(); ++j) {
            new_board[i+upper_x][j] = left_ghost[i][j];
        }
    }

    // middle
    for (int i = 0; i < local_sizex; ++i) {
        for (int j = 0; j < local_sizey; ++j) {
            new_board[i+upper_x][j+left_y] = board[i][j];
        }
    }

    // right
    for(int i = 0; i < right_ghost.size(); ++i) {
        for (int j = 0; j < right_ghost[i].size(); ++j) {
            new_board[i+upper_x][j+left_y+local_sizey] = right_ghost[i][j];
        }
    }

    for (int i = 0; i < lower_left_ghost.size(); ++i) {
        for (int j = 0; j < lower_left_ghost[i].size(); ++j) {
            new_board[i+upper_x+local_sizex][j] = lower_left_ghost[i][j];
        }
    }

    for (int i = 0; i < lower_ghost.size(); ++i) {
        for (int j = 0; j < lower_ghost[i].size(); ++j) {
            new_board[i+upper_x+local_sizex][j+left_y] = lower_ghost[i][j];
        }
    }

    for (int i = 0; i < lower_right_ghost.size(); ++i) {
        for (int j = 0; j < lower_right_ghost[i].size(); ++j) {
            new_board[i+upper_x+local_sizex][j+left_y+local_sizey] = lower_right_ghost[i][j];
        }
    }

    // std::cout << "before calculate all: rank " << rank << std::endl;

    // for (int i = 0; i < new_size; ++i) {
    //     for (int j = 0; j < local_sizey; ++j){
    //         std::cout << new_board[i][j] << " "; 
    //     }
    //     std::cout << std::endl;
    // }

    for (int t = 1 ; t <= update_frequency; ++t) {
        vector<vector<int>> temp_new_board = new_board;

        int x_0 = min(t, upper_x);
        int x_1 = new_x - min(t, lower_x);
        int y_0 = min(t, left_y);
        int y_1 = new_y - min(t, right_y);

        for (int i = x_0; i < x_1; ++i) {
            int row_start = i - 1 >= 0 ? i - 1 :  0;
            int row_end = i + 1 < new_x? i + 1 : new_x - 1;
            for (int j = y_0 ; j < y_1; ++j) {
                int col_start = j - 1 >= 0 ? j - 1 : 0;
                int col_end = j + 1 < new_y ? j + 1 : new_y - 1;
                int alive_neighbour = 0; 
                //inlcude itself
                for (int row = row_start; row <= row_end; row++) {
                    for (int col = col_start; col <= col_end; col++) {
                        alive_neighbour = alive_neighbour + new_board[row][col];
                    }
                }
                alive_neighbour = alive_neighbour - new_board[i][j];
            
                //logistic 
                //self is 1
                if (new_board[i][j] == 1) {
                    if (alive_neighbour == 2 || alive_neighbour == 3){
                    } else {
                        temp_new_board[i][j] = 0;
                    }
                } else {
                    if (alive_neighbour == 3) {
                        temp_new_board[i][j] = 1;
                    }
                }
            } 
        }

        // std::cout << "before calculate all: rank " << rank << std::endl;

        // for (int i = 0; i < new_size; ++i) {
        //     for (int j = 0; j < local_sizey; ++j){
        //         std::cout << temp_new_board[i][j] << " "; 
        //     }
        //     std::cout << std::endl;
        // } 

        new_board = temp_new_board;
    }

    // copy back board values
    for (int i = 0; i < local_sizex; ++i) {
        for (int j = 0; j < local_sizey; ++j) {
            board[i][j] = new_board[i+upper_x][j+left_y];
        }
    }

    // std::cout << "after calculate all: rank " << rank << std::endl;

    // for (int i = 0; i < local_sizex; ++i) {
    //     for (int j = 0; j < local_sizey; ++j){
    //         std::cout << board[i][j] << " "; 
    //     }
    //     std::cout << std::endl;
    // }
    
}


void update(int rank, int step){
    /*send ghost area
    **recevie ghost area
    **decide which rows need to update and calculate 
    **calculate based on time
    */
    vector<vector<int>> upper_ghost;
    vector<vector<int>> lower_ghost;
    vector<vector<int>> left_ghost;
    vector<vector<int>> right_ghost;
    vector<vector<int>> upper_left_ghost;
    vector<vector<int>> upper_right_ghost;
    vector<vector<int>> lower_left_ghost;
    vector<vector<int>> lower_right_ghost;

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

    //for confirm the accuracy of the code, get all ghost area
    // if (step == 1){
    //     std::cout << "rank:" << rank << " upper_ghost:" <<std::endl;
    //     for (int i = 0; i < upper_ghost.size(); i++) {
    //         for (int j = 0; j < local_sizey; ++j){
    //                 std::cout << upper_ghost[i][j] << " "; 
    //         }
    //         std::cout << std::endl;
    //     }

    //     std::cout << "rank:" << rank << " lower_ghost:" <<std::endl;
    //     for (int i = 0; i < lower_ghost.size(); i++) {
    //         for (int j = 0; j < local_sizey; ++j){
    //                 std::cout << lower_ghost[i][j] << " "; 
    //         }
    //         std::cout << std::endl;
    //     }

    // }
    // std::cout << "rank:" << rank << " upper_ghost:" <<std::endl;
    // for (int i = 0; i < upper_ghost.size(); i++) {
    //     for (int j = 0; j < local_sizey; ++j){
    //             std::cout << upper_ghost[i][j] << " "; 
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "rank:" << rank << " lower_ghost:" <<std::endl;
    //  for (int i = 0; i < lower_ghost.size(); i++) {
    //     for (int j = 0; j < local_sizey; ++j){
    //             std::cout << lower_ghost[i][j] << " "; 
    //     }
    //     std::cout << std::endl;
    // }
    

    //TODO:get back
    // calculate_all(rank, step, upper_ghost, lower_ghost);
}






void gather(int rank, int* data, int* disp, int* resvcounts){
    int* recv_temp;
    int* send_data = new int[local_sizex * local_sizey];
    for (int i = 0; i < local_sizex; i++) {
        for (int j = 0; j < local_sizey; j++) {
            send_data[i * local_sizey + j] = board[i][j];
        }
    }
    MPI_Gatherv(send_data, local_sizex * local_sizey, MPI_INT, data, resvcounts, disp, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] send_data;

}