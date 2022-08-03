#!/bin/bash
#SBATCH -J exceptrial
#SBATCH -c 4
#SBATCH --output exceptrial.log
#SBATCH --error exceptrial.err
#SBATCH --mem 40G
#SBATCH --export=ALL 

cd /home/lalonso/projects/goldStd/exception/LSHPGQ
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Bacteroides_fragilis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Enterococcus_faecalis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Escherichia_coli.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Finegoldia_magna.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Fusobacterium_nucleatum.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Gemella_haemolysans.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Helicobacter_pylori.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Lactococcus_lactis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Leptotrichia_goodfellowii.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Listeria_ivanovii.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Mobiluncus_curtisii.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Moraxella_catarrhalis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Mycoplasma_genitalium.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Pasteurella_dagmatis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Propionibacterium_acnes.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Rothia_mucilaginosa.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Staphylococcus_epidermidis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Streptobacillus_moniliformis.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Streptococcus_thermophilus.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Ureaplasma_urealyticum.bam
/home/lalonso/projects/goldStd/scripts/02-BWA-MEM-process_BAM.sh Veillonella_parvula.bam
# 
# LPBQ	done
# HPBQ	done
# HPGQ	done
# LSHPGQ	