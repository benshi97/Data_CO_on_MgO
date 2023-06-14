#!/bin/bash

set -e

datarsync='rsync -zarv --include="*/" --include="*.xyz" --include="*input" --include="*.hist" --include="out" --include="*.output" --include="*EDISP*" --include="OUTCAR" --include="MINP" --include="ED4DISP" --include="INCAR" --include="KPOINTS" --include="POSCAR" --include="CONTCAR" --include="mrcc.out" --include="*.sh" --include="orca.inp" --include="orca.out" --include="orca.bq"  --include="*energ*" --exclude="*"'

for i in 01_PBE-D2-Ne  02_revPBE-D4  03_vdW-DF  04_rev-vdW-DF2  05_PBE0-D4 06_B3LYP-D2-Ne
do
	mkdir -p ${i}
	cd ${i}
	rsync -zarv --include="*/" --include="*.xyz" --include="*input" --include="*.hist" --include="out" --include="*.output" --include="*EDISP*" --include="OUTCAR" --include="MINP" --include="ED4DISP" --include="INCAR" --include="KPOINTS" --include="POSCAR" --include="CONTCAR" --include="mrcc.out" --include="*.sh" --include="orca.inp" --include="orca.out" --include="orca.bq"  --include="*energ*" --exclude="*" archer:/mnt/lustre/a2fs-work3/work/e89/e89/bxs21/Mol-Adsorb/22_09_23-XC_Relax_Compare/02-CO_Surface/${i}/3/ AD_SLAB
	rsync -zarv --include="*/" --include="*.xyz" --include="*input" --include="*.hist" --include="out" --include="*.output" --include="*EDISP*" --include="OUTCAR" --include="MINP" --include="ED4DISP" --include="INCAR" --include="KPOINTS" --include="POSCAR" --include="CONTCAR" --include="mrcc.out" --include="*.sh" --include="orca.inp" --include="orca.out" --include="orca.bq"  --include="*energ*" --exclude="*" archer:/mnt/lustre/a2fs-work3/work/e89/e89/bxs21/Mol-Adsorb/22_09_23-XC_Relax_Compare/01-Surface_Relax/${i}/ SLAB
	rsync -zarv --include="*/" --include="*.xyz" --include="*input" --include="*.hist" --include="out" --include="*.output" --include="*EDISP*" --include="OUTCAR" --include="MINP" --include="ED4DISP" --include="INCAR" --include="KPOINTS" --include="POSCAR" --include="CONTCAR" --include="mrcc.out" --include="*.sh" --include="orca.inp" --include="orca.out" --include="orca.bq"  --include="*energ*" --exclude="*" archer:/mnt/lustre/a2fs-work3/work/e89/e89/bxs21/Mol-Adsorb/22_09_23-XC_Relax_Compare/01-CO_Relax/${i}/ AD


	cd ../
done
