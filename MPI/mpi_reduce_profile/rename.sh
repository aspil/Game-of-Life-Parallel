#!/bin/bash


for f in *.mpiP; do
	n_proc=`echo $f | sed -r 's/.*.x.([0-9]*)\..*/\1/'`
	line=`sed '2q;d' $f`
	mode=`echo $f | sed -r 's/game_([a-z_]*)\.x..*/\1/'`
	grid_size=$(echo $line | tr -dc '0-9')
	# echo $grid_size
	newname=`echo "${mode}_${n_proc}_${grid_size}.mpiP"`
	echo "$f ---> $newname"
	mv $f $newname
done