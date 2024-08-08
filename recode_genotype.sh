#!/bin/bash
#SBATCH --job-name=recode_genotype
#SBATCH --output=recode_genotype_%A_%a.log
#SBATCH --error=recode_genotype_%A_%a.err
#SBATCH --array=1-22   # Adjust this range if you have different chromosome numbers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=standard-cpu

# Ensure the PATH includes the directory where PLINK is installed
export PLINK_DIR=/path/to/plink
export PATH=$PLINK_DIR:$PATH

# Change to the directory where the PLINK files are located
export PLINK_BIN_DIR=/path/to/plink/files
cd $PLINK_BIN_DIR

# Define the chromosome based on the array task ID
CHR=$SLURM_ARRAY_TASK_ID

# Recode genotype data to 0/1/2 format for the specific chromosome without filtering any SNP
plink --bfile chr${CHR} --recode A --out chr${CHR}_recoded --const-fid