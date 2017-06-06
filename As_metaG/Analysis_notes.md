# Analysis of Cen13 metagenomes
### Taylor Dunivin

### Goals
* Check quality and filter Cen13 metagenomic datasets (2 total)
* Test for the presence of AsRGs using Xander (gene targeted metagenome assembler)
* Examine the overlap between gene presence and As resistant isolate collection. 

### Table of contents
* [June 6, 2017]()


### June 6, 2017
Check the quality of Cen13 fastq files: Do this for both fastqc Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.fastq.gz and fastqc Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.fastq.gz
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
