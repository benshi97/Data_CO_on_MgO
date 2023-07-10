#!/bin/bash

set -e

for i in 06_B3LYP-D2-Ne # 01_PBE-D2-Ne 02_revPBE-D4 03_vdW-DF 04_rev-vdW-DF2 05_PBE0-D4 06_B3LYP-D2-Ne
do
	for j in 2.51 2.52 2.53 2.54 2.55 2.56 2.57 2.58 2.60 2.62 # 2.36 2.38 2.40 2.41 2.42 2.43 2.44 2.45 2.46 2.47 2.48 2.49 2.50 2.51 2.52 2.53 2.54 2.55 2.56 2.57 2.58 2.60 2.62 #2.51 2.52 2.53 2.54 2.55 2.56 2.57 2.58 2.60 2.62 #2.36 2.38 2.40 2.41 2.42 2.43 2.44 2.45 2.46 2.47 2.48 2.49 2.50 2.51 2.52 2.53 2.54 2.55 2.56 2.57 2.58 2.60 2.62 
	do
		mkdir -p ${i}/${j}
		cd ${i}/${j}
		cp ../../{KPOINTS,POTCAR} .
		cp ../../INCAR_${i} INCAR

		cp ../../POSCAR_${j} POSCAR
		sed -i 's/IBRION=2/IBRION=-1/g' INCAR
		sed -i 's/NSW=1000/NSW=0/g' INCAR
		cp ../../archer_vasp.sh .
		sbatch archer_vasp.sh
		cd ../../
	done
done
