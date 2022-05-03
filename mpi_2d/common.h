#include <cstdint>
#include <mpi.h>

void initiate(int rank, int sizex,int sizey, int* data, int ranks, int frequency, int rank_row, int rank_col);
void gather(int rank, int sizex, int sizey, int* data, int* disp, int* resvcounts);
void update(int rank,int step, int proc_per_row);