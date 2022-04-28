#include <cstdint>
#include <mpi.h>

void initiate(int rank, int sizex,int sizey, int* data, int ranks);
void gather(int rank, int sizex, int sizey, int* data);