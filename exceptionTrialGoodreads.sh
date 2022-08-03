#!/bin/bash
#SBATCH -J exceptrial
#SBATCH -c 8
#SBATCH --output exceptrial.log
#SBATCH --error exceptrial.err
#SBATCH --mem 40G
#SBATCH --export=ALL 

cd /home/lalonso/projects/goldStd/goodreads/LPBQ
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Bacteroides_fragilis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Enterococcus_faecalis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Escherichia_coli_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Finegoldia_magna_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Fusobacterium_nucleatum_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Gemella_haemolysans_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Helicobacter_pylori_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Lactococcus_lactis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Leptotrichia_goodfellowii_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Listeria_ivanovii_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Mobiluncus_curtisii_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Moraxella_catarrhalis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Mycoplasma_genitalium_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Pasteurella_dagmatis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Cutibacterium_acnes_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Rothia_mucilaginosa_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Staphylococcus_epidermidis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Streptobacillus_moniliformis_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Streptococcus_thermophilus_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Ureaplasma_urealyticum_goodreads_r1.fq.gz
/home/lalonso/projects/goldStd/scripts/02-BWA-aln-process_fastqPE.sh Veillonella_parvula_goodreads_r1.fq.gz
# 
# LPBQ  done
# HPBQ  done
# HPGQ  DONE
# LSHPGQ        

mv *stats ./LSHPGQ
mv *bam* ./LSHPGQ
mv *.tsv ./LSHPGQ


