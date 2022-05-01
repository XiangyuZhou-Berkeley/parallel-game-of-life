#include <cstdint>
#include <mpi.h>

void initiate(int rank, int sizex,int sizey, int* data, int ranks, int frequency = 1);
void gather(int rank, int sizex, int sizey, int* data, int* disp, int* resvcounts);
void update(int rank,int step);