#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ALIVE 'x'
#define DEAD '_'


// void evolve(char** cur, char** next, int pos_i, int pos_j){       // should this return char type? (DEAD/ALIVE)
// 	// dunno yet in what form the grid 'next' will be
// 	int alive_neighbors = 0;
// 	char ALIVE = 'x', DEAD = '-';
	 
// 	if(cells[pos_i][pos_j] == ALIVE){
// 		// 0-1 or 4-8 => die (underpop or overpop)
// 		// 2-3 => live
// 		if(cells[pos_i-1][pos_j-1] == ALIVE) alive_neighbors++;
// 		if(cells[pos_i-1][pos_j] == ALIVE) alive_neighbors++;
// 		if(cells[pos_i-1][pos_j+1] == ALIVE) alive_neighbors++;
// 		if(cells[pos_i][pos_j-1] == ALIVE){
// 			alive_neighbors++;
// 			if(alive_neighbors == 4){
// 				next[pos_i][pos_j] = DEAD;  // die due to overpopulation
// 				return;
// 			}
// 		}
// 		if(cells[pos_i][pos_j+1] == ALIVE){
// 			alive_neighbors++;
// 			if(alive_neighbors >= 4){
// 				next[pos_i][pos_j] = DEAD;  // die due to overpopulation
// 				return;
// 			}
// 		}
// 		if(cells[pos_i+1][pos_j-1] == ALIVE){
// 			alive_neighbors++;
// 			if(alive_neighbors >= 4){
// 				next[pos_i][pos_j] = DEAD;  // die due to overpopulation
// 				return;
// 			}
// 		}
// 		if(cells[pos_i+1][pos_j] == ALIVE){
// 			alive_neighbors++;
// 			if(alive_neighbors >= 4){
// 				next[pos_i][pos_j] = DEAD;  // die due to overpopulation
// 				return;
// 			}
// 		}
// 		if(cells[pos_i+1][pos_j+1] == ALIVE){
// 			alive_neighbors++;
// 			if(alive_neighbors >= 4){
// 				next[pos_i][pos_j] = DEAD;  // die due to overpopulation
// 				return;
// 			}
// 		}

// 		if(alive_neighbors < 2){
// 			next[pos_i][pos_j] = DEAD; // die due to underpopulation
// 			return;
// 		} else {
// 			next[pos_i][pos_j] = ALIVE; // survive => 2 or 3 alive neighbors
// 			return;
// 		}
// 	}
// }


int main(int argc, char *argv[]){
	srand(time(0));

	MPI_Comm cartesian;     // will be the new communicator
	int grid_width = 64, comm_sz;       //maybe grid_width should be passed by arg, e.g 64
	int my_rank, neighbour_rank, src, dest;
	int comm_dims[2];       // dimensions of communicator shape, initialized to 0,0
	int period[2] = {1,1};       // in this implementation, we work on a simple-grid-shaped communicator
	
	int width_local, height_local;      // width x length of this process's chunk, initialized to 0x0
	int chunk_pos[2];       // position of the chunk in the global grid of cells, initialized to 0,0

	MPI_Init(NULL, NULL);   
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Dims_create(comm_sz, 2, comm_dims);    // determines a suitable distribution of processes, according to comm_sz
	MPI_Cart_create(MPI_COMM_WORLD, 2 , comm_dims, period, 0, &cartesian);      // create a new communicator with topology information


	height_local = (int) (grid_width / comm_dims[0]);      // casting is not actually necessary, everything is a power of 2
	width_local = (int) (grid_width / comm_dims[1]);
	
	char **current = malloc((height_local+2) * sizeof(char*));
	char **next =  malloc((height_local+2) * sizeof(char*));

	current[0] = malloc((height_local+2) * (width_local+2) * sizeof(char));
	next[0] = malloc((height_local+2) * (width_local+2) * sizeof(char));

	for (int i = 0; i < (height_local+2); i++){
		current[i] = current[0] + i * (width_local+2);
		next[i] = current[0] + i * (width_local+2);
	}

	// initialize my current subgrid randomly
	// printf("My rank is %d\n", my_rank);
	for (int i = 1; i < height_local + 1; i++) {      // initialize cells; let's see if grid[i][j] works  
		for (int j = 1; j < width_local + 1; j++) {
			current[i][j] = (rand() % 10) < 7 ? DEAD : ALIVE;       // dead : alive ratio 7:3
			// printf("%c", current[i][j]);
		}
		// printf("\n");
	}
	/* Get the coordinates for the current process */
	int my_coords[2];
	MPI_Cart_coords(cartesian, my_rank, 2, my_coords);
	
	/* Precalculate the neighbours' ranks of the current process */
	int neighbours[8];
	int neighbour_dims[2];
	
	MPI_Cart_shift(cartesian, 0, 1, &src, &dest);
	neighbours[0] = src;		/* North neighbour */
	neighbours[1] = dest; 	/* South neighbour */

	MPI_Cart_shift(cartesian, 1, 1, &src, &dest);
	neighbours[2] = src;	/* West neighbour */
	neighbours[3] = dest;	/* East neighbour */

	neighbour_dims[0] = my_coords[0]-1; neighbour_dims[1] = my_coords[1]-1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[4] = neighbour_rank;		/* North-West neighbour */
	
	neighbour_dims[0] = my_coords[0]-1; neighbour_dims[1] = my_coords[1]+1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[5] = neighbour_rank;		/* North-East neighbour */
	
	neighbour_dims[0] = my_coords[0]+1; neighbour_dims[1] = my_coords[1]-1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[6] = neighbour_rank;		/* South-West neighbour */
	
	neighbour_dims[0] = my_coords[0]+1; neighbour_dims[1] = my_coords[1]+1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[7] = neighbour_rank;		/* South-East neighbour */

	MPI_Datatype column;		//we construct a new type to send the columns
	MPI_Type_vector(height_local, 1, width_local + 2, MPI_CHAR, &column); // count, blocklength, stride, old, new
	MPI_Type_commit(&column);
	
	MPI_Request SRequest[8], RRequest[8];
	/* Northern process communication */
	MPI_Send_init(&current[1][1], width_local, MPI_CHAR, neighbours[0], 1, cartesian, &SRequest[0]);
	MPI_Recv_init(&current[0][1], 			1, MPI_CHAR, neighbours[0], 1, cartesian, &RRequest[0]);
	
	/* Southern process communication */
	MPI_Send_init(&current[height_local][1],   width_local, MPI_CHAR, neighbours[1], 2, cartesian, &SRequest[1]);
	MPI_Recv_init(&current[height_local+1][1], 			 1, MPI_CHAR, neighbours[1], 2, cartesian, &RRequest[1]);

	/* Western process communication */
	MPI_Send_init(&current[1][1], height_local, column, neighbours[2], 3, cartesian, &SRequest[2]);
	MPI_Recv_init(&current[1][0], 			 1, column, neighbours[2], 3, cartesian, &RRequest[2]);
	
	/* Easter process communication */
	MPI_Send_init(&current[1][width_local],  height_local, column, neighbours[3], 4, cartesian, &SRequest[3]);
	MPI_Recv_init(&current[1][width_local+1], 			1, column, neighbours[3], 4, cartesian, &RRequest[3]);

	/* North-western process communication */
	MPI_Send_init(&current[1][1], 1, MPI_CHAR, neighbours[4], 5, cartesian, &SRequest[4]);
	MPI_Recv_init(&current[0][0], 1, MPI_CHAR, neighbours[4], 5, cartesian, &RRequest[4]);
	
	/* North-eastern process communication */
	MPI_Send_init(&current[width_local][1],	  1, MPI_CHAR, neighbours[5], 6, cartesian, &SRequest[5]);
	MPI_Recv_init(&current[0][width_local+1], 1, MPI_CHAR, neighbours[5], 6, cartesian, &RRequest[5]);
	
	/* South-western process communication */
	MPI_Send_init(&current[height_local][1],   1, MPI_CHAR, neighbours[6], 7, cartesian, &SRequest[6]);
	MPI_Recv_init(&current[height_local+1][0], 1, MPI_CHAR, neighbours[6], 7, cartesian, &RRequest[6]);
	
	/* South-eastern process communication */
	MPI_Send_init(&current[height_local][width_local], 	   1, MPI_CHAR, neighbours[7], 8, cartesian, &SRequest[7]);
	MPI_Recv_init(&current[height_local+1][width_local+1], 1, MPI_CHAR, neighbours[7], 8, cartesian, &RRequest[7]);
	

	// MPI_Comm_free(&cartesian);
	printf("my rank = %d\n",my_rank);
	MPI_Finalize();
	printf("Now running serial\n");
	
	free(current[0]);
	free(current);
	free(next[0]);
	free(next);
	return 0;
}