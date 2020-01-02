#!/bin/bash
#SBATCH --job-name=bbmap_simple  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sawagz01@nyu.edu # Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=80gb # Job memory request
#SBATCH --time=09:00:00 # Time limit hrs:min:sec
#SBATCH --array=1,2,3
#SBATCH --output=/gpfs/scratch/sawagz01/Bioinformatics/Assignment3/kallisto_%j.log # Standard output and error log
#SBATCH -p cpu_short

module load kallisto

kallisto quant -i transcripts_HG38.idx -o BU${SLURM_ARRAY_TASK_ID}_kallisto -b 100 /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU${SLURM_ARRAY_TASK_ID}_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlBU${SLURM_ARRAY_TASK_ID}_2.fastq.gz

kallisto quant -i transcripts_HG38.idx -o CtrlDs12H${SLURM_ARRAY_TASK_ID}_kallisto -b 100 /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H${SLURM_ARRAY_TASK_ID}_1.fastq.gz /gpfs/data/courses/bmscga2604/Assignment1/CtrlDs12H${SLURM_ARRAY_TASK_ID}_2.fastq.gz
