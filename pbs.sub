#!/bin/bash
#SBATCH --job-name=JD15_6
#SBATCH -A aosc-hi
#SBATCH -t 12:00:00
##SBATCH -N 2
#SBATCH -n 16
#SBATCH --share
#SBATCH --mail-user=tchor@umd.edu
#SBATCH --mail-type=ALL
#SBATCH --constraint="rhel6"

module purge
module load hdf5/intel/2016.3.210/intelmpi/shared/1.10.0
ulimit -s unlimited

cd `pwd`

/bin/echo "mpirun -np $SLURM_NTASKS diablo_dev/diablo starting at `date`"
mpirun -np $SLURM_NTASKS diablo_dev/diablo > output${SLURM_JOBID}.dat
