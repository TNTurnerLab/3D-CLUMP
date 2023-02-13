#!/bin/bash
#SBATCH --job-name=GI_GW  # Job name
#SBATCH --ntasks=16
#SBATCH --mem=48gb           # Job memory request
#SBATCH --time=72:00:00        # Time limit hrs:min:sec
#SBATCH --output=/mnt/disk005/data/projects/GRUMP/log/GI_GWjob.log
echo "test"

source /home/ychen/.bashrc
source activate snakemake
snakemake -c8  -s /mnt/disk005/data/projects/GRUMP/CLUMP-master/high_throughput_GI_GW/case.control.snake --unlock
snakemake -c8  -s /mnt/disk005/data/projects/GRUMP/CLUMP-master/high_throughput_GI_GW/case.control.snake
