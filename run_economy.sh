#!/bin/bash
### Project code
#PBS -A UMCP0012
### Job Name
#PBS -N JD15_2d
### Merge output and error files
#PBS -j oe
#PBS -q economy

### Select X nodes with Y CPUs and Y mpiprocs each for a total of X*Y MPI processes
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=12:00:00
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M tomaschor@gmail.com
#PBS -o log.out

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR
export MPI_USE_ARRAY=false

### Run the executable
#module load impi
source /glade/u/home/tomasc/opt/load_modules_diablo.sh
mpirun -np 16 source/diablo > output.log

