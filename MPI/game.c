#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define ALIVE 1
#define DEAD 0
#define MAX_GENERATIONS 100


int main(int argc, char *argv[]){
	srand(time(0));

	MPI_Comm cartesian;     // will be the new communicator
	int grid_size, comm_sz;       //maybe grid_size should be passed by arg, e.g 64
	int my_rank, neighbour_rank, src, dest;
	int comm_dims[2];       // dimensions of communicator shape, initialized to 0,0
	int period[2] = {1,1};       // in this implementation, we work on a simple-grid-shaped communicator
	
	int width_local, height_local;      // width x length of this process's chunk, initialized to 0x0
	int i, j;       // declaring here, so it only costs once
	// int chunk_pos[2];       // position of the chunk in the global grid of cells, initialized to 0,0


#ifdef REDUCE
    int similarity_counter = 0, empty_counter = 0, sim_sum, e_sum;
#endif
	if (argc != 2) {
		printf("Wrong program arguments. Usage: mpirun [-n] [X] game.x <grid size>\n");
		return EXIT_FAILURE;
	}
	grid_size = atoi(argv[1]);

	MPI_Init(NULL, NULL);   
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Dims_create(comm_sz, 2, comm_dims);    // determines a suitable distribution of processes, according to comm_sz
	MPI_Cart_create(MPI_COMM_WORLD, 2 , comm_dims, period, 0, &cartesian);      // create a new communicator with topology information

	height_local = (int) (grid_size / comm_dims[0]);      // casting is not actually necessary, everything is a power of 2
	width_local = (int) (grid_size / comm_dims[1]);
	/* Allocate the 'before' and 'after' grids */
	// printf("EDW1\n");
	char **current = malloc((height_local+2) * sizeof(char*));
	char *current0 = malloc((height_local+2) * (width_local+2) * sizeof(char));

	char **next =  malloc((height_local+2) * sizeof(char*));
	char *next0 = malloc((height_local+2) * (width_local+2) * sizeof(char));

	for (i = 0; i < (height_local+2); i++){
		current[i] = current0 + i * (width_local+2);
		next[i] = next0 + i * (width_local+2);
	}

	/* Initialize the 'before' or current subgrid */
	for (i = 1; i < height_local + 1; i++)
		for (j = 1; j < width_local + 1; j++)
			current[i][j] = (rand() % 10) < 8 ? DEAD : ALIVE;	// dead : alive ratio 7:3

	/* Get the coordinates for the current process */
	int my_coords[2];
	MPI_Cart_coords(cartesian, my_rank, 2, my_coords);
	
	/* Precalculate the neighbours' ranks of the current process */
	int neighbours[8];
	int neighbour_dims[2];
	
	MPI_Cart_shift(cartesian, 0, 1, &src, &dest);
	neighbours[0] = src;	/* North neighbour */
	neighbours[1] = dest; 	/* South neighbour */

	MPI_Cart_shift(cartesian, 1, 1, &src, &dest);
	neighbours[2] = src;	/* West neighbour */
	neighbours[3] = dest;	/* East neighbour */
	
	/* North-West neighbour */
	neighbour_dims[0] = my_coords[0]-1; neighbour_dims[1] = my_coords[1]-1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[4] = neighbour_rank;		
	
	/* North-East neighbour */
	neighbour_dims[0] = my_coords[0]-1; neighbour_dims[1] = my_coords[1]+1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[5] = neighbour_rank;		
	
	/* South-West neighbour */
	neighbour_dims[0] = my_coords[0]+1; neighbour_dims[1] = my_coords[1]-1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[6] = neighbour_rank;		
	
	/* South-East neighbour */
	neighbour_dims[0] = my_coords[0]+1; neighbour_dims[1] = my_coords[1]+1;
	MPI_Cart_rank(cartesian, neighbour_dims, &neighbour_rank);
	neighbours[7] = neighbour_rank;		

	MPI_Datatype column;		//we construct a new type to send the columns
	MPI_Type_vector(height_local, 1, width_local + 2, MPI_CHAR, &column); // count, blocklength, stride, old, new
	MPI_Type_commit(&column);
	
	MPI_Request SRequests[8];
	MPI_Request RRequests[8];
	/* Northern process communication */
	MPI_Recv_init(&(current[0][1]), width_local, MPI_CHAR, neighbours[0], 0, cartesian, &RRequests[0]);
	MPI_Send_init(&(current[1][1]), width_local, MPI_CHAR, neighbours[0], 0, cartesian, &SRequests[0]);

	/* Southern process communication */
	MPI_Recv_init(&(current[height_local+1][1]), width_local, MPI_CHAR, neighbours[1], 0, cartesian, &RRequests[1]);
	MPI_Send_init(&(current[height_local][1]),   width_local, MPI_CHAR, neighbours[1], 0, cartesian, &SRequests[1]);

	/* Western process communication */
	MPI_Recv_init(&(current[1][0]), 1, column, neighbours[2], 0, cartesian, &RRequests[2]);
	MPI_Send_init(&(current[1][1]), 1, column, neighbours[2], 0, cartesian, &SRequests[2]);
	
	/* Easter process communication */
	MPI_Recv_init(&(current[1][width_local+1]), 1, column, neighbours[3], 0, cartesian, &RRequests[3]);
	MPI_Send_init(&(current[1][width_local]),  	1, column, neighbours[3], 0, cartesian, &SRequests[3]);

	/* North-western process communication */
	MPI_Recv_init(&(current[0][0]), 1, MPI_CHAR, neighbours[4], 0, cartesian, &RRequests[4]);
	MPI_Send_init(&(current[1][1]), 1, MPI_CHAR, neighbours[4], 0, cartesian, &SRequests[4]);
	
	/* North-eastern process communication */
	MPI_Recv_init(&(current[0][width_local+1]), 1, MPI_CHAR, neighbours[5], 0, cartesian, &RRequests[5]);
	MPI_Send_init(&(current[1][width_local]),   1, MPI_CHAR, neighbours[5], 0, cartesian, &SRequests[5]);
	
	/* South-western process communication */
	MPI_Recv_init(&(current[height_local+1][0]), 1, MPI_CHAR, neighbours[6], 0, cartesian, &RRequests[6]);
	MPI_Send_init(&(current[height_local][1]),   1, MPI_CHAR, neighbours[6], 0, cartesian, &SRequests[6]);
	
	/* South-eastern process communication */
	MPI_Recv_init(&(current[height_local+1][width_local+1]), 1, MPI_CHAR, neighbours[7], 0, cartesian, &RRequests[7]);
	MPI_Send_init(&(current[height_local][width_local]), 	 1, MPI_CHAR, neighbours[7], 0, cartesian, &SRequests[7]);
	
	char **temp;
	int alive_neighbors = 0;
	int generation = 0;
	MPI_Barrier(cartesian);
	double time_start = MPI_Wtime();
	MPI_Pcontrol(1);

	while (generation < MAX_GENERATIONS) {
		MPI_Startall(8, RRequests);
		MPI_Startall(8, SRequests);
		/* Εvolution of inner (white) cells */
		for(i = 2; i < height_local; i++){
			for(j = 2; j < width_local; j++){
				alive_neighbors = current[i-1][j-1] +
								  current[i-1][j+1] +
								  current[i][j-1]	+
								  current[i-1][j]	+
								  current[i][j+1]   +
								  current[i+1][j-1] +
								  current[i+1][j]   +
								  current[i+1][j+1];
				if(current[i][j] == ALIVE)
					next[i][j] = ((alive_neighbors == 2) || (alive_neighbors == 3)) ? ALIVE : DEAD;
				else
					next[i][j] = (alive_neighbors == 3) ? ALIVE : DEAD;
#ifdef REDUCE
				similarity_counter += current[i][j] ^ next[i][j];	// if 0 at the end of loop, then similar
				empty_counter += next[i][j];		// if 0 then none is alive => empty
#endif
			}
		}

		MPI_Waitall(8,  RRequests, MPI_STATUSES_IGNORE);

		for(i = 1; i < height_local + 1; i++){   // calculating the green columns
			/*-------------  Leftmost column -------------*/
			alive_neighbors = current[i-1][0] +
							  current[i-1][1] +
							  current[i-1][2] +
							  current[i][0]   +
							  current[i][2]   +
							  current[i+1][0] +
							  current[i+1][1] +
							  current[i+1][2];

			if(current[i][1] == ALIVE)
				next[i][1] = ((alive_neighbors == 2) || (alive_neighbors == 3)) ? ALIVE : DEAD;
			else
				next[i][1] = (alive_neighbors == 3) ? ALIVE : DEAD;
#ifdef REDUCE
			similarity_counter += current[i][1] ^ next[i][1];	// if 0 at the end of loop, then similar
			empty_counter += next[i][1];		// if 0 then none is alive => empty
#endif
			/*-------------  Rightmost column -------------*/
			alive_neighbors = current[i-1][width_local-1] +
							  current[i-1][width_local]   +
							  current[i-1][width_local+1] +
							  current[i][width_local-1]   +
							  current[i][width_local+1]   +
							  current[i+1][width_local-1] +
							  current[i+1][width_local]   +
							  current[i+1][width_local+1];

			if(current[i][width_local] == ALIVE)
				next[i][width_local] = ((alive_neighbors == 2) || (alive_neighbors == 3)) ? ALIVE : DEAD;
			else
				next[i][width_local] = (alive_neighbors == 3) ? ALIVE : DEAD;
#ifdef REDUCE
			similarity_counter += current[i][width_local] ^ next[i][width_local];	// if not 0 at the end of loop, then similar
			empty_counter += next[i][width_local];		// if 0 then none is alive => empty
#endif
		}

		for(j = 1; j < width_local + 1; j++){   // calculating the green rows
			/* -------------- Uppermost row --------------*/
			alive_neighbors = current[0][j-1] +
							  current[0][j]   +
							  current[0][j+1] +
							  current[1][j-1] +
							  current[1][j+1] +
							  current[2][j-1] +
							  current[2][j]   +
							  current[2][j+1];

			if(current[1][j] == ALIVE)
				next[1][j] = ((alive_neighbors == 2) || (alive_neighbors == 3)) ? ALIVE : DEAD;
			else
				next[1][j] = (alive_neighbors == 3) ? ALIVE : DEAD;
#ifdef REDUCE
			similarity_counter += current[1][j] ^ next[1][j];	// if 0 at the end of loop, then similar
			empty_counter += next[1][j];		// if 0 then none is alive => empty
#endif
			/* -------------- Lowermost row --------------*/
			alive_neighbors = current[height_local-1][j-1] +
							  current[height_local-1][j]   +
							  current[height_local-1][j+1] +
							  current[height_local][j-1]   +
							  current[height_local][j+1]   +
							  current[height_local+1][j-1] +
							  current[height_local+1][j]   +
							  current[height_local+1][j+1];

			if(current[height_local][j] == ALIVE)
				next[height_local][j] = ((alive_neighbors == 2) || (alive_neighbors == 3)) ? ALIVE : DEAD;
			else
				next[height_local][j] = (alive_neighbors == 3) ? ALIVE : DEAD;
#ifdef REDUCE
			similarity_counter += current[height_local][j] ^ next[height_local][j];	// if 0 at the end of loop, then similar
			empty_counter += next[height_local][j];		// if 0 then none is alive => empty
#endif
		}

		MPI_Waitall(8,  SRequests, MPI_STATUSES_IGNORE);
		generation++; 
		temp = current;
		current = next;
		next = temp;

#ifdef REDUCE
		if (generation % 10 == 0) {
			MPI_Allreduce(&similarity_counter, &sim_sum, 1, MPI_INT, MPI_SUM, cartesian);
			MPI_Allreduce(&empty_counter, &e_sum , 1, MPI_INT, MPI_SUM, cartesian);
			// if (sim_sum == 0) {
			// 	// Similar grids
			// }
			// if (!e_sum) {
			// 	// Empty grid
			// }
		}
#endif      
	}	/* End of main loop */
	MPI_Pcontrol(0);
	double time_end = MPI_Wtime();
	if(my_rank == 0) {
		printf("Time elapsed = %f\n", time_end- time_start);
	}
	MPI_Type_free(&column);
	MPI_Comm_free(&cartesian);
	MPI_Finalize();
	
	// free(current[0]);
	// free(current);
	// free(next[0]);
	// free(next);
	return 0;
}