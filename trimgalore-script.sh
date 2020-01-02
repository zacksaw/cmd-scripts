#!/bin/bash
#SBATCH --job-name=bbmap_simple  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=# Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=09:00:00 # Time limit hrs:min:sec
#SBATCH --array=1,2,3
#SBATCH --output=/gpfs/scratch/sawagz01/Bioinformatics/Assignment2/bbmap_simple_{$SLURM_ARRAY_TASK_ID}_%j.log # Standard output and error log
#SBATCH -p cpu_short

### SECTION #2 - TRIM DATA ###

module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here!!!!
module load fastqc/0.11.7

trim_galore --paired --length 30 -o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/ /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H${SLURM_ARRAY_TASK_ID}_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H${SLURM_ARRAY_TASK_ID}_2.fastq.gz --fastqc --fastqc_args "-o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/" --q 30 --gzip

trim_galore --paired --length 30 -o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/ /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU${SLURM_ARRAY_TASK_ID}_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU${SLURM_ARRAY_TASK_ID}_2.fastq.gz --fastqc --fastqc_args "-o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/" --q 30 --gzip


