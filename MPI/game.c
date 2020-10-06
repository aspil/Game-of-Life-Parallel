#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


char** create_subgrid(const int dims[]){
    srand(time(0));
    char **subgrid = malloc(dims[0] * sizeof(char));      //allocating space for grid
    char ALIVE = 'x', DEAD = '-';

    for (int i = 0; i < dims[0]; i++){
        subgrid[i] = malloc(dims[1] * sizeof(char));
        for(int j = 0; j < dims[1]; j++){
            subgrid[i][j] = (rand() % 10) < 7 ? DEAD : ALIVE;       // dead : alive ratio 7:3
        }
    }
    return subgrid;
}


void evolve(char** cur** next, int pos_i, int pos_j){       // should this return char type? (DEAD/ALIVE)
    // dunno yet in what form the grid 'next' will be
    int alive_neighbors = 0;
    char ALIVE = 'x', DEAD = '-';
     
    if(cells[pos_i][pos_j] == ALIVE){
        // 0-1 or 4-8 => die (underpop or overpop)
        // 2-3 => live
        if(cells[pos_i-1][pos_j-1] == ALIVE) alive_neighbors++;
        if(cells[pos_i-1][pos_j] == ALIVE) alive_neighbors++;
        if(cells[pos_i-1][pos_j+1] == ALIVE) alive_neighbors++;
        if(cells[pos_i][pos_j-1] == ALIVE){
            alive_neighbors++;
            if(alive_neighbors == 4){
                next[pos_i][pos_j] = DEAD;  // die due to overpopulation
                return;
            }
        }
        if(cells[pos_i][pos_j+1] == ALIVE){
            alive_neighbors++;
            if(alive_neighbors >= 4){
                next[pos_i][pos_j] = DEAD;  // die due to overpopulation
                return;
            }
        }
        if(cells[pos_i+1][pos_j-1] == ALIVE){
            alive_neighbors++;
            if(alive_neighbors >= 4){
                next[pos_i][pos_j] = DEAD;  // die due to overpopulation
                return;
            }
        }
        if(cells[pos_i+1][pos_j] == ALIVE){
            alive_neighbors++;
            if(alive_neighbors >= 4){
                next[pos_i][pos_j] = DEAD;  // die due to overpopulation
                return;
            }
        }
        if(cells[pos_i+1][pos_j+1] == ALIVE){
            alive_neighbors++;
            if(alive_neighbors >= 4){
                next[pos_i][pos_j] = DEAD;  // die due to overpopulation
                return;
            }
        }

        if(alive_neighbors < 2){
            next[pos_i][pos_j] = DEAD; // die due to underpopulation
            return;
        } else {
            next[pos_i][pos_j] = ALIVE; // survive => 2 or 3 alive neighbors
            return;
        }  
}


int main(void){
    int grid_width = 64, comm_sz;    //maybe grid_width should be passed by arg, e.g 64
    int my_rank;
    int chunk_dims[2] = {0,0};      // width x length of this process's chunk, initialized to 0x0
    int chunk_pos[2] = {0,0};     // position of the chunk in the global grid of cells, initialized to 0,0
    int comm_dims[2] = {0,0};       // dimensions of communicator shape, initialized to 0,0
    int period[2] = {0,0};       // in this implementation, we work on a simple-grid-shaped communicator
    int my_coords[2];       // position of this process in the communicator grid

    MPI_Comm cartesian;     // will be the new communicator

    MPI_Init(NULL, NULL);   
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    MPI_Dims_create(comm_sz, 2, comm_dims);    // determines a suitable distribution of processes, according to comm_sz
    MPI_Cart_create(MPI_COMM_WORLD, 2 , comm_dims, period, 0, &cartesian);      // create a new communicator with topology information

    chunk_dims[0] = (int) (grid_width / comm_dims[0]) + 4;      // casting is not actually necessary, everything is a power of 2
    chunk_dims[1] = (int) (grid_width / comm_dims[1]) + 4;

    ////////////// beginning of distinguished process behavior //////////////

    MPI_Cart_Coords(cartesian, my_rank, 2, my_coords);      // retrieve this process's position in the cartesian comm
    
    // calculation of the chunk position for this process
    // chunk_pos[1] = my_coords[1] * chunk_dims[1];
    // chunk_pos[0] = my_coords[0] * chunk_dims[0];

    // create my part of the grid, locally
    char **my_grid = create_subgrid(chunk_dims);
    
    // do we need to construct the type of 'row' / 'column' ??
    // switch(rank){
        // MPI_IRecv to neighbors
        // MPI_ISend to neighbors
    //}
    


}