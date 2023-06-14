#!/bin/bash

set -e

for i in 2.51 2.53 2.55 2.56 2.57 2.58 2.60 2.62 # 2.36 2.38 2.40 2.41 2.42 2.43 2.44 2.45 2.46 2.47 2.48 2.49 2.50 2.52 2.54
do
	for j in MP2 CC
	do
		for k in TZ QZ
		do
			for l in AD_SLAB SLAB_CP AD_CP
			do
				mkdir -p ${i}/${j}/${k}/${l}
				cd ${i}/${j}/${k}/${l}
				cp ../../../../Input/DF${j}_CO_rdf_EMBD_bond_dist_${i}_${k}_MINP_${l} MINP
				cp ../../../../cirrus_mrcc.sh .
				sbatch cirrus_mrcc.sh
				cd ../../../../
			done
		done
	done
done

