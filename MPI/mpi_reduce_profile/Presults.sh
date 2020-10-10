#!/bin/bash

for f in *.mpiP; do
	# Get the number of processes
	n_proc=`echo $f | sed -r 's/.*_([0-9]*)_.*.mpiP/\1/'`

	# Get the execution type
	mode=`echo  $f | sed -r 's/([a-z]*)_.*/\1/'`
	# Get the size of the grid
	line=`sed '2q;d'  $f`
	grid_size=$(echo $line | tr -dc '0-9')

	echo "Mode: $mode, Processes: $n_proc, Grid Size: $grid_size"
	offset=$((3+$n_proc))
	cat $f | grep -A $offset -e "MPI Time (seconds)" | tail -n 1
	echo
done