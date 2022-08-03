#!/bin/bash
### input is BAM from TCGA non human reads to keep the cluster jobs easy ### in lava
# ssh lava 
cd /scratch/lalonso/exception/HPGQ

# 0- sort reads by name in the bam not necessary because with all the preprocessing finally they become unsorted
# 1- convert bam to fastq
# 2- trim adapters or remove reads with adapters with skewer.
# 3- remove low complexity and mean low qual   

for inputbam in `ls ./*.bam`; do 
  prefixinput=`echo $inputbam | sed 's/.bam//i' | sed 's/\.\///i' `
  goodreads=`echo $prefixinput | sed 's/$/_goodreads/i' `
  samtools fastq $inputbam --threads 8 -N 2>error_samtools.tmp | skewer - --min 48 -n --quiet -1 --threads 8 --mode pe -x /opt/skewer-master/adapters/adapters.fa 2>error_skewer.tmp | perl /opt/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -min_len 48 -out_format 3 -out_good stdout -out_bad null -no_qual_header -min_qual_mean 10 -ns_max_p 5 -lc_method dust -lc_threshold 10 -trim_tail_left 5 -trim_tail_right 5 -trim_ns_left 1 -trim_ns_right 1 -trim_qual_left 10 -trim_qual_right 10 2> $goodreads.prinseq.log | gzip > $goodreads.fq.gz
done

for gzipped in `ls *_goodreads.fq.gz`; do 
  /scratch/lalonso/exception/single_paired_filter_fj.py $gzipped
done

