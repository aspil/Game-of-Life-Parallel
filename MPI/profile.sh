#!/bin/bash
for((sz=320;sz<=320*2;sz=2*sz)); do
	for((x=1;x<=64;x=2*x)); do
		echo "Running with $x processes on $sz x $sz grid"
		./sub.sh $x $sz
	done
done