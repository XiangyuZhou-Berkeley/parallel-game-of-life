#include <mpi.h>
#include <iostream>
#include <vector>
using namespace std;


int local_sizex = 0;
int local_sizey = 0;
int final_sizex = 0;
int update_frequency = 1;
int total_rank = 0;

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

void rearrange(vector<vector<int>> &after, int* before, int sizex, int sizey){
    // cout << "vector size   " << after.size() << std::endl;
    if (after.size() != 0) {
        int temp_count = 0;
        for (int i = 0; i < sizex; i++) {
            vector<int> temp_list;
            for (int j = 0; j < sizey; j++) {
                temp_list.push_back(before[temp_count++]);
            }
            after.push_back(temp_list);
        }
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


void update(int rank, int step, int proc_per_row){
    /*send ghost area
    **recevie ghost area
    **decide which rows need to update and calculate 
    **calculate based on time
    */
    int proc_row = rank / proc_per_row;
    int proc_col = rank % proc_per_row;


    vector<vector<int>> upper_ghost;
    vector<vector<int>> lower_ghost;
    vector<vector<int>> left_ghost;
    vector<vector<int>> right_ghost;
    vector<vector<int>> upper_left_ghost;
    vector<vector<int>> upper_right_ghost;
    vector<vector<int>> lower_left_ghost;
    vector<vector<int>> lower_right_ghost;

    int num_transfer_row = update_frequency * local_sizey;
    int num_transfer_corner = update_frequency * update_frequency;
    int num_transfer_col = update_frequency * local_sizey;

    int* s_upper_ghost = (int*) malloc((int)num_transfer_row  * sizeof(int));
    int* s_lower_ghost = (int*) malloc((int)num_transfer_row  * sizeof(int));
    int* r_upper_ghost = (int*) malloc((int)num_transfer_row  * sizeof(int));
    int* r_lower_ghost = (int*) malloc((int)num_transfer_row  * sizeof(int));


    int* s_left_ghost = (int*) malloc((int)num_transfer_col  * sizeof(int));
    int* s_right_ghost = (int*) malloc((int)num_transfer_col  * sizeof(int));
    int* r_left_ghost = (int*) malloc((int)num_transfer_col  * sizeof(int));
    int* r_right_ghost = (int*) malloc((int)num_transfer_col  * sizeof(int));

    int* s_upper_left_ghost = (int*) malloc((int)num_transfer_corner * sizeof(int));
    int* s_upper_right_ghost = (int*) malloc((int)num_transfer_corner  * sizeof(int));
    int* s_lower_left_ghost = (int*) malloc((int)num_transfer_corner * sizeof(int));
    int* s_lower_right_ghost = (int*) malloc((int)num_transfer_corner  * sizeof(int));
    int* r_upper_left_ghost = (int*) malloc((int)num_transfer_corner * sizeof(int));
    int* r_upper_right_ghost = (int*) malloc((int)num_transfer_corner  * sizeof(int));
    int* r_lower_left_ghost = (int*) malloc((int)num_transfer_corner * sizeof(int));
    int* r_lower_right_ghost = (int*) malloc((int)num_transfer_corner  * sizeof(int));






    MPI_Request requests[8]; 
    // 0 1 2 
    // 3   4
    // 5 6 7  
    MPI_Status statuses[8];
    int num_send_row = min(update_frequency, local_sizex);
    int num_send_col = min(update_frequency, local_sizey);
 
    int num_received_row_ghost = 0;
    int num_received_col_ghost = 0;
    int num_received_corner_ghost = 0;

    int num_send_row_ghost = num_send_row * local_sizey;
    int num_send_col_ghost = num_send_col * local_sizex;
    int num_send_corner_ghost = update_frequency * update_frequency;

    //send ghost area
    
    //only s_lower_ghost
    if (proc_row == 0) {
        int num_temp = 0;
        for (int i = local_sizex - num_send_row; i < local_sizex; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_lower_ghost[num_temp++] = board[i][j];
            }
        }
        MPI_Isend(s_lower_ghost, num_send_row_ghost, MPI_INT, rank+proc_per_row, rank, MPI_COMM_WORLD, &requests[6]);
    } else if (proc_row == proc_per_row-1) {
        //only s_upper_ghost,possible here only a few
        int num_temp = 0;
        for (int i = 0; i < num_send_row; i++) {
            for (int j = 0; j < local_sizey; j++) {
                s_upper_ghost[num_temp++] = board[i][j];
            }
        }
        // send upper ghost to rank (rank - proc_per_row)
        MPI_Isend(s_upper_ghost, num_send_row_ghost, MPI_INT, rank-proc_per_row, rank, MPI_COMM_WORLD, &requests[1]);
    } else {
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
        MPI_Isend(s_upper_ghost, num_send_row_ghost, MPI_INT, rank-proc_per_row, rank, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(s_lower_ghost, num_send_row_ghost, MPI_INT, rank+proc_per_row, rank, MPI_COMM_WORLD, &requests[6]);    
    }

    /*
    send left and right ghost 
    */
    if (proc_col == 0) {
        // only send right
        int num_temp = 0;
        for (int i = 0; i < local_sizex; ++i ) {
            for (int j = local_sizey - num_send_col; j < local_sizey; ++j) {
                s_right_ghost[num_temp++] = board[i][j];
            }
        }
        // send right ghost to proc ( rank + 1)
        MPI_Isend(s_right_ghost, num_send_col_ghost, MPI_INT, rank+1, rank, MPI_COMM_WORLD, &requests[4]);    
    } else if (proc_col == proc_per_row - 1) {
        // only send left
        int num_temp = 0;
        for (int i = 0; i < local_sizex; ++i ) {
            for (int j = 0; j < num_send_col; ++j) {
                s_left_ghost[num_temp++] = board[i][j];
            }
        } 
        // send left ghost to proc (rank - 1)
        MPI_Isend(s_left_ghost, num_send_col_ghost, MPI_INT, rank-1, rank, MPI_COMM_WORLD, &requests[3]);  
    } else {
        // send right
        int num_temp = 0;
        for (int i = 0; i < local_sizex; ++i ) {
            for (int j = local_sizey - num_send_col; j < local_sizey; ++j) {
                s_right_ghost[num_temp++] = board[i][j];
            }
        }

        num_temp = 0;
        for (int i = 0; i < local_sizex; ++i ) {
            for (int j = 0; j < num_send_col; ++j) {
                s_left_ghost[num_temp++] = board[i][j];
            }
        } 
        MPI_Isend(s_left_ghost, num_send_col_ghost, MPI_INT, rank-1, rank, MPI_COMM_WORLD, &requests[3]);  
        MPI_Isend(s_right_ghost, num_send_col_ghost, MPI_INT, rank+1, rank, MPI_COMM_WORLD, &requests[4]);   
    }

    /* 
    send corners ghost
    */
    if (proc_col == 0 ){
        if (proc_row == 0) {
            // send lower right 
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_lower_right_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_lower_right_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row+1, rank, MPI_COMM_WORLD, &requests[7]);  
        } else if (proc_row == proc_per_row - 1 ){
            // send upper right
            int num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_upper_right_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_upper_right_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row+1, rank, MPI_COMM_WORLD, &requests[2]);  
        } else {
            // send lower right and upper right
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_lower_right_ghost[num_temp++] = board[i][j];
                }
            }
            num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_upper_right_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_upper_right_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row+1, rank, MPI_COMM_WORLD, &requests[2]); 
            MPI_Isend(s_lower_right_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row+1, rank, MPI_COMM_WORLD, &requests[7]);  
        }
    } else if (proc_col == proc_per_row - 1 ) {
        if (proc_row == 0) {
            // send lower left
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_lower_left_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_lower_left_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row-1, rank, MPI_COMM_WORLD, &requests[5]); 
        } else if (proc_row == proc_per_row - 1 ){
            // send upper left
            int num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_upper_left_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_upper_left_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row-1, rank, MPI_COMM_WORLD, &requests[0]); 
        } else {
            // send lower left and upper left
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_lower_left_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_upper_left_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_upper_left_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row-1, rank, MPI_COMM_WORLD, &requests[0]); 
            MPI_Isend(s_lower_left_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row-1, rank, MPI_COMM_WORLD, &requests[5]); 
        } // else
    } else {
         if (proc_row == 0) {
            // send lower left and lower right
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_lower_left_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_lower_right_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_lower_left_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row-1, rank, MPI_COMM_WORLD, &requests[5]); 
            MPI_Isend(s_lower_right_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row+1, rank, MPI_COMM_WORLD, &requests[7]); 

        } else if (proc_row == proc_per_row - 1 ){
            // send upper left and upper right
            int num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_upper_left_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_upper_right_ghost[num_temp++] = board[i][j];
                }
            }
            MPI_Isend(s_upper_left_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row-1, rank, MPI_COMM_WORLD, &requests[0]); 
            MPI_Isend(s_upper_right_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row+1, rank, MPI_COMM_WORLD, &requests[2]); 

        } else {
            // send all
            int num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_lower_left_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = local_sizex - update_frequency; i < local_sizex; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_lower_right_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = 0; j < update_frequency; ++j) {
                    s_upper_left_ghost[num_temp++] = board[i][j];
                }
            }

            num_temp = 0;
            for (int i = 0; i < update_frequency; ++i) {
                for (int j = local_sizey - update_frequency; j < local_sizey; ++j) {
                    s_upper_right_ghost[num_temp++] = board[i][j];
                }
            }

            MPI_Isend(s_upper_left_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row-1, rank, MPI_COMM_WORLD, &requests[0]); 
            MPI_Isend(s_upper_right_ghost, num_send_corner_ghost, MPI_INT, rank-proc_per_row+1, rank, MPI_COMM_WORLD, &requests[2]); 
            MPI_Isend(s_lower_left_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row-1, rank, MPI_COMM_WORLD, &requests[5]); 
            MPI_Isend(s_lower_right_ghost, num_send_corner_ghost, MPI_INT, rank+proc_per_row+1, rank, MPI_COMM_WORLD, &requests[7]); 

        } // else 
    }

    /* receive begin
       0 1 2
       3   4
       5 6 7
    */

    if (proc_row == 0) {

        //receive lower
        MPI_Recv(r_lower_ghost,num_send_row_ghost, MPI_INT, rank + proc_per_row, rank + proc_per_row, MPI_COMM_WORLD, &statuses[6]);
        lower_ghost.resize(num_send_row_ghost);
        if (proc_col == 0) {

            //receive right and lower_right
            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_lower_right_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row + 1, rank + proc_per_row + 1, MPI_COMM_WORLD, &statuses[7]);

            right_ghost.resize(num_send_col_ghost);
            lower_right_ghost.resize(num_send_corner_ghost);
        } else if (proc_col == proc_per_row - 1) {

            //receive left and receive lower_left
            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_lower_left_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row - 1, rank + proc_per_row - 1, MPI_COMM_WORLD, &statuses[5]);

            left_ghost.resize(num_send_col_ghost);
            lower_left_ghost.resize(num_send_corner_ghost);
        } else {

            //receive right and lower_right  left and receive lower_left
            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_lower_right_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row + 1, rank + proc_per_row + 1, MPI_COMM_WORLD, &statuses[7]);
            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_lower_left_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row - 1, rank + proc_per_row - 1, MPI_COMM_WORLD, &statuses[5]);

            right_ghost.resize(num_send_col_ghost);
            lower_right_ghost.resize(num_send_corner_ghost);
            left_ghost.resize(num_send_col_ghost);
            lower_left_ghost.resize(num_send_corner_ghost);
        }
    } else if (proc_row == proc_per_row - 1) {

        //receive upper
        MPI_Recv(r_upper_ghost,num_send_row_ghost, MPI_INT, rank - proc_per_row, rank - proc_per_row, MPI_COMM_WORLD, &statuses[1]);

        upper_ghost.resize(num_send_row_ghost);
        
        if (proc_col == 0) {

            //receive right and upper_right
            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_upper_right_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row + 1, rank - proc_per_row + 1, MPI_COMM_WORLD, &statuses[2]);
            
            right_ghost.resize(num_send_col_ghost);
            upper_right_ghost.resize(num_send_corner_ghost);
        } else if (proc_col == proc_per_row - 1) {

            //receive left and receive upper_left
            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_upper_left_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row - 1, rank - proc_per_row - 1, MPI_COMM_WORLD, &statuses[0]);

            left_ghost.resize(num_send_col_ghost);
            upper_left_ghost.resize(num_send_corner_ghost);
        } else {

            //receive right and upper_right     left and receive upper_left
            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_upper_right_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row + 1, rank - proc_per_row + 1, MPI_COMM_WORLD, &statuses[2]);
            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_upper_left_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row - 1, rank - proc_per_row - 1, MPI_COMM_WORLD, &statuses[0]);

            right_ghost.resize(num_send_col_ghost);
            upper_right_ghost.resize(num_send_corner_ghost);
            left_ghost.resize(num_send_col_ghost);
            upper_left_ghost.resize(num_send_corner_ghost);
        }

    } else {

        //recv upper    lower
        MPI_Recv(r_upper_ghost,num_send_row_ghost, MPI_INT, rank - proc_per_row, rank - proc_per_row, MPI_COMM_WORLD, &statuses[1]);
        MPI_Recv(r_lower_ghost,num_send_row_ghost, MPI_INT, rank + proc_per_row, rank + proc_per_row, MPI_COMM_WORLD, &statuses[6]);

        upper_ghost.resize(num_send_row_ghost);
        lower_ghost.resize(num_send_row_ghost);
        
        if (proc_col == 0) {
            //recv right, upper_right, lower_right
            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_upper_right_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row + 1, rank - proc_per_row + 1, MPI_COMM_WORLD, &statuses[2]);
            MPI_Recv(r_lower_right_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row + 1, rank + proc_per_row + 1, MPI_COMM_WORLD, &statuses[7]);

            right_ghost.resize(num_send_col_ghost);
            upper_right_ghost.resize(num_send_corner_ghost);
            lower_right_ghost.resize(num_send_corner_ghost);
        } else if (proc_col == proc_per_row - 1) {

            //recv left, upper_left, lower_left
            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_upper_left_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row - 1, rank - proc_per_row - 1, MPI_COMM_WORLD, &statuses[0]);
            MPI_Recv(r_lower_left_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row - 1, rank + proc_per_row - 1, MPI_COMM_WORLD, &statuses[5]);

            left_ghost.resize(num_send_col_ghost);
            upper_left_ghost.resize(num_send_corner_ghost);
            lower_left_ghost.resize(num_send_corner_ghost);
        } else {

            //recv right, upper_right, lower_right  left, upper_left, lower_left

            MPI_Recv(r_right_ghost,num_send_col_ghost, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD, &statuses[4]);
            MPI_Recv(r_upper_right_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row + 1, rank - proc_per_row + 1, MPI_COMM_WORLD, &statuses[2]);
            MPI_Recv(r_lower_right_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row + 1, rank + proc_per_row + 1, MPI_COMM_WORLD, &statuses[7]);

            MPI_Recv(r_left_ghost,num_send_col_ghost, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &statuses[3]);
            MPI_Recv(r_upper_left_ghost,num_send_corner_ghost, MPI_INT, rank - proc_per_row - 1, rank - proc_per_row - 1, MPI_COMM_WORLD, &statuses[0]);
            MPI_Recv(r_lower_left_ghost,num_send_corner_ghost, MPI_INT, rank + proc_per_row - 1, rank + proc_per_row - 1, MPI_COMM_WORLD, &statuses[5]);

            right_ghost.resize(num_send_col_ghost);
            upper_right_ghost.resize(num_send_corner_ghost);
            lower_right_ghost.resize(num_send_corner_ghost);
            left_ghost.resize(num_send_col_ghost);
            upper_left_ghost.resize(num_send_corner_ghost);
            lower_left_ghost.resize(num_send_corner_ghost);
        }

    }

    

    // if (rank == 0) {
    //     for (int i = 0; i < local_sizey * update_frequency;i++) {
    //        cout << r_lower_ghost[i] << " ";
    //     }
    //     cout << std::endl;
    // }
    
    //put into grid
    rearrange(upper_left_ghost, r_upper_left_ghost, update_frequency, update_frequency);
    rearrange(upper_ghost, r_upper_ghost, update_frequency, local_sizey);
    rearrange(upper_right_ghost, r_upper_right_ghost, update_frequency, update_frequency);
    rearrange(left_ghost, r_left_ghost, local_sizex, update_frequency);
    rearrange(right_ghost, r_right_ghost, local_sizex, update_frequency);
    rearrange(lower_left_ghost, r_lower_left_ghost, update_frequency, update_frequency);
    rearrange(lower_ghost, r_lower_ghost, update_frequency, local_sizey);
    rearrange(lower_right_ghost, r_lower_right_ghost, update_frequency, update_frequency);

    
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
    calculate_all(rank, step, upper_ghost, lower_ghost, left_ghost, right_ghost, upper_left_ghost, upper_right_ghost, lower_left_ghost, lower_right_ghost);
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