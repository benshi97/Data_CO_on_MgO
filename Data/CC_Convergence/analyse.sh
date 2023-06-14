#!/bin/bash

set -e

touch time_sheet_cc
rm time_sheet_cc

for i in 1_Canonical  1_Local  2_Local  3_Local
do
	for j in DZ TZ
	do
		for k in AD_SLAB SLAB_CP AD_CP
		do
			cd ${i}/${j}/${k}
			a=$(head -n 20 mrcc.out | tail -n 1 | awk ' { print $2 " " $3 } ')
			b=$(tail -n 3 mrcc.out | head -n 1 | awk ' { print $2 " " $3 } ')
			echo "${a}  ${b}" >> ../../../time_sheet_cc
			cd ../../../
		done
	done
done
