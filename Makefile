exec = mpi_test.x

MPI/exec:
	mpicc MPI/mpi_test.c -o MPI/mpi_test.x

.PHONY:
	clean

clean:
	rm *.o* *.e*