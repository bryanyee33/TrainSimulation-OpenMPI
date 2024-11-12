#!/bin/bash

## You can use sbatch instead of salloc to run batch jobs, as usual.

#SBATCH --job-name=a3-mpi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=xs-4114
#SBATCH --time=00:03:00
#SBATCH --output=a3_%j.slurmlog
#SBATCH --error=a3_%j.slurmlog

echo "Received mpirun arguments: $@"

# Display commands being run in stdout
set -x

echo "Seq:"
mpirun --report-bindings --display-map --display-allocation --map-by node --bind-to core \
./bench_seq testcases/performance/min_test.in

echo "Done"