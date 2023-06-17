#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --time=72:00:00
## Number of MPI tasks (MPI)
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=40

#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=2000MB
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --partition=l_long
#SBATCH --qos=ll



##Set OMP_NUM_THREADS to the same value as -c
##with a fallback in case it isn't set.
##SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
 omp_threads=$SLURM_CPUS_PER_TASK
else
 omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads
echo "OMP_NUM_THREADS=" $omp_threads

#========================================
# load modules and run simulation

module unload PrgEnv-cray
module load PrgEnv-intel

license=5-1563
releasename=544_18Apr17
variant=""
build=cnl7.0_intel19.1.3.304
releasever=5.4.4


bindir="/ddn/projects/vasp/5-2728/544_18Apr17/cnl6.0_intel17.0.1.132/src/vasp.5.4.4/bin"
srun --propagate=STACK,MEMLOCK --hint=nomultithread  $bindir/vasp_std
