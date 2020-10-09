#!/bin/bash

n_proc=$1
grid_sz=$2
directive=$3
if [ $n_proc -le 0 ]; then
	echo "[$(whoami)]: Please provide a positive number of processes."
	exit
fi

res=$(echo "l($n_proc)/l(2)" | bc -l)
res_ceil=$(echo "$res" | perl -nl -MPOSIX -e 'print ceil($_);')
res_floor=$(echo "$res" | perl -nl -MPOSIX -e 'print floor($_);')

if [[ $res_ceil != $res_floor ]]; then
	echo "[$(whoami)]: $1 is not a power of 2. Please give a number that is a power of 2."
	exit
fi

# Remove the switch argument from mpirun command on mpiPBSscript.sh
sed -i -e 's/mpirun -n [0-9]* game[a-zA-Z_]*.x [0-9]*/mpirun game_'$directive'.x '$2'/g' mpiPBSscript.sh

if [[ n_proc -le 4 ]]; then
	sed -i -e 's/select=[0-9]*/select=1/g' mpiPBSscript.sh
	sed -i -e 's/mpiprocs=[0-9]*/mpiprocs='$n_proc'/g' mpiPBSscript.sh
	sed -i -e 's/mpirun game[a-zA-Z_]*.x [0-9]*/mpirun -n '$n_proc' game_'$directive'.x '$grid_sz'/g' mpiPBSscript.sh
else
	nodes=$(($n_proc / 8))
	sed -i -e 's/select=[0-9]*/select='$nodes'/g' mpiPBSscript.sh
	sed -i -e 's/mpiprocs=[0-9]*/mpiprocs=8/g' mpiPBSscript.sh
fi
sed -i -e 's/#PBS -N myJob_[0-9]*_[0-9]*/#PBS -N myJob_'$n_proc'_'$grid_sz'/g' mpiPBSscript.sh
cat mpiPBSscript.sh | grep mpirun
echo
# Finally sumbit the job on the PBS queue
# qsub mpiPBSscript.sh