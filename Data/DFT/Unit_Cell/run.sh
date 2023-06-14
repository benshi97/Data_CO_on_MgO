#!/bin/bash

set -e

for i in 01_PBE-D2-Ne  02_revPBE-D4  03_vdW-DF  04_rev-vdW-DF2  05_PBE0-D4  06_B3LYP-D2-Ne
	do
		cd ${i}
		cp ../archer_vasp.sh .
		sbatch archer_vasp.sh
		cd ../
	done

