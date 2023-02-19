#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --time=168:00:00
#SBATCH --partition=long7
#SBATCH --ntasks=40
#SBATCH --exclusive
#SBATCH --hint=nomultithread
##SBATCH --mem-per-cpu=2000



#=========================================
# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads
echo "OMP_NUM_THREADS=" $omp_threads

#========================================
# load modules and run simulation
#export LAMMPS_USE=_omp
#module load lammps/17Nov16
module unload PrgEnv-intel/6.0.9
module load  PrgEnv-cray/6.0.9 
module load vasp/vasp5/5.4
ulimit -s unlimited

bindir="/lustre/software/vasp/vasp5/vasp.5.4.4.pl2/bin/"
srun   $bindir/vasp_std
