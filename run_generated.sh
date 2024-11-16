#!/bin/bash

## You can use sbatch instead of salloc to run batch jobs, as usual.

#SBATCH --job-name=a3-mpi
#SBATCH --nodes=1
#SBATCH --ntasks=2
##SBATCH --partition=xs-4114
#SBATCH --partition=i7-13700
#SBATCH --time=00:03:00
#SBATCH --output=gen/a3_%j.slurmlog
#SBATCH --error=gen/a3_%j.slurmlog

echo "We are running on $(hostname)"
echo "Received mpirun arguments: $@"

# Display commands being run in stdout
set -x

touch output.txt
> output.txt # empty file
cd testcases

for folder in line_length num_ticks num_trains stations
do  
    cd $folder
    {
    for infile in $(ls)
    do
        mpirun --map-by node --bind-to core ../../trains $infile >> ../../output.txt
        echo
    done
    }
    echo
    cd ..
done

echo "Done"