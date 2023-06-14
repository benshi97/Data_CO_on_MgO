#!/bin/bash

set -e

for i in 1 2 3
do
	for j in 1 2
	do
		mkdir -p ${i}x${i}x1/${j}
		cd ${i}x${i}x1/${j}
		cp ../../../01-PWSCF/${j}/pw.d .
		cp ../../../01-PWSCF/${j}/*.upf .
		sed -i 's/K_POINTS GAMMA/K_POINTS automatic/g' pw.d
		echo -e "\n${i} ${i} 1 0 0 0" >> pw.d
		cp ../../csd3_qe.sh .
		sbatch csd3_qe.sh
		cd ../../
	done
done
