#!/bin/bash
#SBATCH -J mgo_relax
#SBATCH -A T2-CS152-CPU
#SBATCH -p icelake
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=38
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00
#SBATCH --mail-type=NONE
##SBATCH --no-requeue

numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
module load intel-oneapi-mkl/2022.1.0/intel/mngj3ad6

application="/home/bxs21/Programs/icelake/vasp.6.3.0/bin/vasp_std"
options=""
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

export OMP_NUM_THREADS=1

np=$[${numnodes}*${mpi_tasks_per_node}]

export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"
$CMD
