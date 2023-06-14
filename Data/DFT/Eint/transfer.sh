#!/bin/bash

set -e

for i in 01_PBE-D2-Ne  02_revPBE-D4  03_vdW-DF  04_rev-vdW-DF2  05_PBE0-D4 06_B3LYP-D2-Ne
do
	mkdir -p ${i}
	cd ${i}
	scp -r diskhost:/zfs/scox_michaelides/michaelides/backup/benshi/2022_CAM_MOL/Final_Data/DFT/04-Int_Ene/CO/${i}/02_revPBE-D4/{AD_SLAB,SLAB_FS,AD_FS} .
	cd ../
done
