#!/bin/bash

## You can use sbatch instead of salloc to run batch jobs, as usual.

#SBATCH --job-name=a3-mpi
#SBATCH --nodes=1
#SBATCH --ntasks=8
##SBATCH --partition=xs-4114
#SBATCH --partition=i7-13700
#SBATCH --time=00:03:00
#SBATCH --output=a3_%j.slurmlog
#SBATCH --error=a3_%j.slurmlog

echo "Received mpirun arguments: $@"

# Display commands being run in stdout
set -x

echo "Mine:"
# Run mpirun with debugging options to display mapping and binding
mpirun --report-bindings --display-map --display-allocation --map-by node --bind-to core \
./trains testcases/performance/hard.in

echo "Done"