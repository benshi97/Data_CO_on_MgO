#!/bin/bash -l
#$ -N vasp
#$ -P Gold
#$ -A UKCP_CAM_C
#$ -l h_rt=48:00:00
#$ -l mem=4700M
#$ -pe mpi 200
#$ -cwd

module unload -f compilers mpi
module load gcc-libs
module load compilers/intel/2019/update5
module load mpi/intel/2019/update5/intel
module load python3/recommended


# 9. Run our MPI job. GERun is a wrapper that launches MPI jobs on Legion.
gerun /home/mmm0606/Programs/vasp.6.3.0/bin/vasp_std >> vasp_output.$JOB_ID
