#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	MPI_Comm comm_cart;
	int rank, comm_sz, id;
	
	int dims[2] = {0,0}, period[2] = {1,1};
	int coords[2];

	MPI_Init(&argc , &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);

	MPI_Dims_create(comm_sz, 2, dims);

	MPI_Cart_create(MPI_COMM_WORLD, 2 , dims, period, 0, &comm_cart);
	printf("dims = %d, %d\n", dims[0], dims[1]);
	if (rank == 0) {
		printf("comm_sz = %d, my_rank = %d\n", comm_sz, rank);
		// coord[0]=0; coord[1]=3;
		// MPI_Cart_rank(comm_cart, coord, &id);
		// printf("The processor at position (%d, %d) has rank %d\n", coord[0], coord[1], id);fflush(stdout);
	}
	else {
		printf("my_rank = %d\n", rank);
	}
		
	MPI_Finalize();
	
	return 0;
}






