#!/bin/bash

# Get parameters from user

read -p 'How many increased volume runs? (+1% volume ea.) ' numrunsinc
read -p 'How many decreased volume runs? (-1% volume ea.) ' numrunsdec
read -p 'Name of structure? ' structname
read -p 'Name of starting folder containing relaxed structure? ' imagefolder 

echo 'Generating increased volume folders...'

for ((i = 1; i <= $numrunsinc; i++)); do
	dir=./$structname-$i-percent-larger
	mkdir $dir
	cp -r ./$imagefolder/* $dir/
	latticeparam=$(sed -n '2p' $dir/POSCAR)
	if (( $i >= 10 )); then
		newlatticeparam=$(echo "print($latticeparam * 1.$i)" | python)
	else
		newlatticeparam=$(echo "print($latticeparam * 1.0$i)" | python)
	fi
	echo $newlatticeparam
	sed -i "2d" $dir/POSCAR
	sed -i "2i $newlatticeparam" $dir/POSCAR
	# also update the INCAR file IBRION tag
	sed -i "3d" $dir/INCAR
	sed -i "3i IBRION=-1" $dir/INCAR
done

echo 'Generating decreased volume folders...'

for ((j = 1; j <= $numrunsdec; j++)); do
	dir=./$structname-$j-percent-smaller
	mkdir $dir
	cp -r ./$imagefolder/* $dir/
	latticeparam=$(sed -n '2p' $dir/POSCAR)
	newlatticeparam=$(echo "print($latticeparam * (1.00 - (0.01 * $j)))" | python)
	echo $newlatticeparam
	sed -i "2d" $dir/POSCAR
	sed -i "2i $newlatticeparam" $dir/POSCAR
	# also update the INCAR file IBRION tag
	sed -i "3d" $dir/INCAR
	sed -i "3i IBRION=-1" $dir/INCAR
done

# generate new kpoints based on updated POSCARs
# also running vasp interactively here since these calculations are so fast...
for d in $structname-*/; do
	# skip the relaxation calculation, we already did it
	if [[ $d == $imagefolder ]]; then
		continue
	fi
	cd $d
	~/bin/kpoints.x
	# this is for running interactively
	#~/bin/vasp6_serial

	# this is for running as seperate jobs
	sbatch runjob.sh
	echo "Job $d submitted"
	cd ..
done

echo "All jobs submitted, wait for completion."
