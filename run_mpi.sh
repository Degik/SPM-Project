#!/bin/bash
#SBATCH --job-name=wavefront_mpi_job
#SBATCH --output=log.out
#SBATCH --error=err.err
#SBATCH --partition=normal
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00

# Set the total number of processes
TOTAL_TASKS=$(( $SLURM_NNODES * $SLURM_NTASKS_PER_NODE ))

echo "Esecuzione del programma MPI su $SLURM_NNODES nodi con $SLURM_NTASKS_PER_NODE task per nodo (totale $TOTAL_TASKS task)."

# Use mpirun on the total tasks
srun --mpi=pmix -n $TOTAL_TASKS ./wavefront_mpi 1000
