#!/bin/bash
#SBATCH --job-name=bbmap_simple  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= # Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=80gb # Job memory request
#SBATCH --time=12:00:00 # Time limit hrs:min:sec
#SBATCH --array=1,2,3
#SBATCH --output=/gpfs/scratch/sawagz01/Bioinformatics/Assignment2/bbmap_simple_{$SLURM_ARRAY_TASK_ID}_%j.log # Standard output and error log
#SBATCH -p cpu_short

module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here!!!!
module load fastqc/0.11.7
module load bbmap/38.25
module load samtools/1.3
module load bedtools/2.26.0

samtools view -b -o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sam
samtools sort -o /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.bam
samtools index /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam

## Extract reads originating from fragments mapping to forward (top) strand (NOTE: THIS IS ONLY APPLICABLE TO dUTP STRANDED RNA-SEQ METHODS)

samtools view -b -f99 /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}fwd2.bam
samtools view -b -f147 /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}fwd1.bam
samtools merge -f /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.forward.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}fwd1.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}fwd2.bam

## Extract reads originating from fragments mapping to reverse (bottom) strand (NOTE: THIS IS ONLY APPLICABLE TO dUTP STRANDED RNA-SEQ METHODS)

samtools view -b -f83 /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}rev2.bam
samtools view -b -f163 /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}rev1.bam
samtools merge -f /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.reverse.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}rev1.bam /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}rev2.bam

## Index files

samtools index /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.reverse.bam
samtools index /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.forward.bam

## Parse using Bedtools (genomeCoverageBed) to generate a simple four-column bedgraph file that can be easily loaded into Gviz

samtools view -b /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.forward.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/data/courses/bmscga2604/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.forward.bedgraph

samtools view -b /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.reverse.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/data/courses/bmscga2604/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.reverse.bedgraph

samtools view -b /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}out.sorted.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/data/courses/bmscga2604/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sawagz01/Bioinformatics/Assignment2/BU${SLURM_ARRAY_TASK_ID}sorted.bedgraph

