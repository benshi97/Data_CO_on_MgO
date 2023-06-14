#!/bin/bash

set -e

for i in 1 2 3 4
do
	for j in AD_SLAB SLAB AD
	do
		mkdir -p ${i}/${j}
		cd ${i}/${j}
		cp ../../Input/${i}/POSCAR_${j} POSCAR
		cp ../../Input/${i}/KPOINTS .
		cp ../../Input/${i}/POTCAR_${j} POTCAR
		cp ../../Input/${i}/INCAR_${j} INCAR
		cp ../../archer_vasp.sh .
		sbatch archer_vasp.sh
		cd ../../
	done
done
