#!/bin/bash
#SBATCH -J mgo_relax
#SBATCH -A T2-CS152-CPU
#SBATCH -p icelake
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --tasks-per-node=19
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
##SBATCH --no-requeue

numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')


export OMP_NUM_THREADS=4

module purge
module load rhel8/default-icl



module load quantum-espresso/7.0/intel/intel-oneapi-mkl/intel-oneapi-mpi/glfi6dek

srun pw.x -i ./pw.d > outpw
