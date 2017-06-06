# Analysis of Cen13 metagenomes
### Taylor Dunivin

### Goals
* Check quality and filter Cen13 metagenomic datasets (2 total)
* Test for the presence of AsRGs using Xander (gene targeted metagenome assembler)
* Examine the overlap between gene presence and As resistant isolate collection. 

### Table of contents
* [June 6, 2017]()


### June 6, 2017
1. Check the quality of Cen13 fastq files: Do this for both fastqc Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.fastq.gz and fastqc Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.fastq.gz
```
#!/bin/bash -login
 
### define resources needed:
### walltime - how long you expect the job to run
#PBS -l walltime=04:00:00
 
### nodes:ppn - how many nodes & cores per node (ppn) that you require
#PBS -l nodes=1:ppn=1
 
### mem: amount of memory that the job will need
#PBS -l mem=100gb
 
### you can give your job a name for easier identification
#PBS -N R1
 
### load necessary modules, e.g.
module load GNU/4.4.5
module load FastQC/0.11.3
 
### change to the working directory where your code is located
cd /mnt/scratch/dunivint/c13
 
### call your executable
fastqc Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.fastq.gz
```

Output: 
[R1](https://github.com/ShadeLab/Arsenic_Growth_Analysis/blob/master/As_metaG/data/Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001_fastqc.html): Quality looks good, but will still trim
[R2](https://github.com/ShadeLab/Arsenic_Growth_Analysis/blob/master/As_metaG/data/Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001_fastqc.html): Quality is not ideal. Will trim

2. Quality trim data files using fastx
```
#load modules
module load GNU/4.4.5
module load FASTX/0.0.14

#quality filter
fastq_quality_filter -Q33 -q 30 -p 50 -z -i Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.fastq -o Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.fastq.gz
fastq_quality_filter -Q33 -q 30 -p 50 -z -i Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.fastq -o Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.fastq.gz
```
