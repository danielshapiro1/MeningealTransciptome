#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 72:00:00
#SBATCH --output=scRNAseq_meninges_alignment_TBI
#SBATCH --partition=standard
#SBATCH -A lukenslab

module load cellranger/3.1.0
cellranger count --id=run_meninges_072820_TBI --transcriptome=/scratch/ds3vz/LukensLab-92447355/FASTQ_Generation_2020-07-21_13_33_44Z-285512227/refdata-gex-mm10-2020-A/ --fastqs=/scratch/ds3vz/LukensLab-92447355/FASTQ_Generation_2020-07-21_13_33_44Z-285512227 --sample=TBI
