#!/bin/bash
#SBATCH --job-name=vcf_to_plink
#SBATCH --output=vcf_to_plink_%A_%a.log
#SBATCH --error=vcf_to_plink_%A_%a.err
#SBATCH --array=1-22   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --partition=standard-cpu

# Ensure the PATH includes the directory where PLINK is installed
export PLINK_DIR=/path/to/plink
export PATH=$PLINK_DIR:$PATH

# Change to the directory where the VCF files are located
export VCF_DIR=/path/to/vcf/files
cd $VCF_DIR

# Define the chromosome based on the array task ID
CHR=$SLURM_ARRAY_TASK_ID

# Convert VCF to PLINK binary format for the specific chromosome
plink --vcf chr${CHR}.dose.filtered.R2_0.8.vcf.gz --make-bed --out chr${CHR} --const-fid
