#!/bin/bash

VASPBIN="path to vasp_gam binary"
CC4SBIN="path to Cc4s binary"

VASP="mpirun ?? $VASPBIN"
CC4S="mpirun ?? $CC4SBIN"

#cc4s.in, POSCAR, POTCAR and KPOINTS need to be present

enc=500
egw=300

cat POSCAR
rm WAVECAR

cat >KPOINTS<<!
Automatically generated mesh
       0
gamma
 1 1 1
!



echo "++++++++++++++++++++++++++++++"
echo "RUN DFT to get a converged guess for HF"
echo "++++++++++++++++++++++++++++++"

cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-6
NCORE=2
!
cat INCAR
$VASP
cp OUTCAR OUTCAR.DFT


echo "++++++++++++++++++++++++++++++"
echo "RUN HF"
echo "++++++++++++++++++++++++++++++"

cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-6
LHFCALC=.TRUE.
AEXX=1.0
ALGO=C
NCORE=2
!
cat INCAR
$VASP
cp OUTCAR OUTCAR.HF


nb=`awk <OUTCAR.HF "/maximum number of plane-waves:/ { print \\$5*2-1 }"`


echo "++++++++++++++++++++++++++++++"
echo "RUN HF diag"
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-6
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO = sub ; NELM = 1
NBANDS = $nb
NCORE=2
!
cat INCAR
$VASP
cp OUTCAR OUTCAR.HFdiag
cp WAVECAR WAVECAR.diag

cp WAVECAR.diag WAVECAR

echo "++++++++++++++++++++++++++++++"
echo "RUN MP2 NOs"
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO = MP2NO ;
NBANDS = $nb
LAPPROX=.TRUE.
!
rm WAVEDER
cat INCAR
$VASP
cp OUTCAR OUTCAR.MP2-NOs


nocc=`awk <OUTCAR.HF "/NELEC/ { print \\$3/2 }"`

echo "going to use $nb bands"
echo "there are $nocc occupied bands"


for nbfp in 101
do


echo "for fpX we need " $nbfp
nbno=`awk <OUTCAR.HF "/NELEC/ { print (\\$3/2)*$nbfp }"`
echo "for fpX we need $nbno  bands"

cp WAVECAR.FNO WAVECAR

echo "++++++++++++++++++++++++++++++"
echo "RUN HF diag of NOs using " $nbno " bands."
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-6
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO = sub ; NELM = 1
NBANDS = $nbno
NBANDSHIGH = $nbno
NCORE=2
!
rm WAVEDER
cat INCAR
$VASP
cp OUTCAR OUTCAR.HFdiag-NOs


echo "++++++++++++++++++++++++++++++"
echo "RUN MP2"
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO = MP2 
NBANDS = $nbno
NBANDSHIGH = $nbno
LSFACTOR=.TRUE.
!
rm WAVEDER
cat INCAR
$VASP
cp OUTCAR OUTCAR.MP2-CBS.$nbfp

done


for nbfp in 11
do


echo "for fpX we need " $nbfp
nbno=`awk <OUTCAR.HF "/NELEC/ { print (\\$3/2)*$nbfp }"`
echo "for fpX we need $nbno  bands"

cp WAVECAR.FNO WAVECAR

echo "++++++++++++++++++++++++++++++"
echo "RUN HF diag of NOs using " $nbno " bands."
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-6
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO = sub ; NELM = 1
NBANDS = $nbno
NBANDSHIGH = $nbno
NCORE=2
!
rm WAVEDER
cat INCAR
$VASP
cp OUTCAR OUTCAR.HFdiag-NOs


echo "++++++++++++++++++++++++++++++"
echo "Dump CC4S input using " $nbno " bands."
echo "++++++++++++++++++++++++++++++"


cat >INCAR <<!
ENCUT = $enc
SIGMA=0.0001
EDIFF = 1E-5
LHFCALC=.TRUE.
AEXX=1.0
ISYM=-1
ALGO=CC4S
NBANDS = $nbno
NBANDSHIGH = $nbno
ENCUTGW=$egw
ENCUTGWSOFT=$egw
!
cat INCAR
$VASP
cp OUTCAR OUTCAR.CC4S

#In this step the following files will be written that are needed for CC4S
#
#
#FockOperator.yaml, FockOperator.dat
#GridVectors.yaml, GridVectors.dat
#CoulombPotential.yaml, CoulombPotential.dat
#DeltaPPHH.yaml, DeltaPPHH.dat
#DeltaHH.yaml, DeltaHH.dat
#CoulombVertexSingularVectors.yaml, CoulombVertexSingularVectors.dat
#CoulombVertex.yaml, CoulombVertex.dat
#EigenEnergies.yaml, EigenEnergies.dat
#Spins.yaml, Spins.dat

echo "++++++++++++++++++++++++++++++"
echo "Run CC4S using " $nb " bands."
echo "++++++++++++++++++++++++++++++"

rm atrip-checkpoint.yaml
$CC4S -i cc4s.in | tee  cc4s.stdout.$nbno.bands


done
