#!/bin/bash
#SBATCH -J exceptrial
#SBATCH -c 4
#SBATCH --output exceptrial.log
#SBATCH --error exceptrial.err
#SBATCH --mem 60G
#SBATCH --export=ALL 

cd /home/lalonso/projects/bladder/bamfiles1
## (1.2) input is BAM, sequencing was Paired End but unmapped reads were singletons ### 
/home/lalonso/projects/bladder/scripts/01-align_bladderUnmapped.sh U0001.unmapped.bam