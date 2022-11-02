#!/bin/bash
#
#SBATCH --job-name=gen_dashboards
#SBATCH --output=logs/%A_%x.log
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=00:30:00


# Get arguments
while getopts n:p: flag
do
    case "${flag}" in
        n) name=${OPTARG};;
        p) password=${OPTARG};;
    esac
done

# Generate dashboards
ml R/4.0
Rscript "scripts/${name}.R"

# Encrypt dashboards
module purge
find "outputs/${name}" -type f -name "*.html" -exec staticrypt {} $password -o {} -r 1 \;

# cleanup
squeue -j $SLURM_JOBID --Format=TimeUsed
