mpi:
	mpicc -Wall -o game.x game.c 

mpi_reduce:
	mpicc -D REDUCE -Wall -o game.x game.c

openmp:
	# mpicc -D openmp -Wall -o game.x game.c

.PHONY:
	clean_mpi, clean_omp

clean_mpi:
	rm -f *.o* *.e* game.x

clean_omp:
