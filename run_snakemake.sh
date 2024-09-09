#!/bin/bash 
#SBATCH --job-name=metagenomics #Job name to be displayed by for example squeue 
#SBATCH --output=metagenomics_%j.out #Path to the file where the job (error) output is written to 
#SBATCH --nodes=1 #Number of nodes. Multiple nodes are only useful for jobs with distributed memory (e.g. MPI). 
#SBATCH --ntasks-per-node=40 #Number of (MPI) processes per node. More than one useful only  for MPI jobs. Maximum number depends nodes (number of cores) 
#SBATCH --mem=128GB  #Memory (RAM) per node. Number followed by unit prefix. 
#SBATCH --partition=smp #Partition/queue in which run the job. 
echo "SLURM_JOBID="$SLURM_JOBID 
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST 
echo "SLURM_NNODES"=$SLURM_NNODES 
echo "SLURMTMPDIR="$SLURMTMPDIR 
echo "Date = $(date)"           
echo "Hostname = $(hostname -s)"        
echo "" 
echo "Number of Nodes Allocated       = $SLURM_JOB_NUM_NODES" 
echo "Number of Tasks Allocated       = $SLURM_NTASKS" 
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK" 
echo "Working Directory = $(pwd)" 
echo "working directory = "$SLURM_SUBMIT_DI 

# Start time
start_time=$(date +%s)

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh  # Adjust the path if Conda is installed elsewhere
conda activate snakemake

# Change directory
cd /lustre/pulkit.h/snakemake_local/gutmicrobiome/pipeline || exit

snakemake --unlock
snakemake --cores all --forceall

echo workflow completed

# End time
end_time=$(date +%s)

# Calculate and print the elapsed time
elapsed_time=$((end_time - start_time))
echo "Total time taken: $elapsed_time seconds"

/bin/hostname |tee result 