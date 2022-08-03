#!/bin/bash
### Remove putative remaining adapters, remove Low Complexity reads, trim Low Quality ends, remove too short reads
# input is a BAM file from TCGA non-human (non-host) reads 
# input may also be a FASTQfile

cd /scratch/lola/exception/HPGQ
# Recursively perform the following actions over a set of BAM files (already mapped to a host genome) in a directory
# 1- convert bam to fastq
# 2- trim adapters or remove reads with adapters with skewer.
# 3- remove low complexity and mean low qual with prinseq
# 4- save good reads in fastq.gz that will contain separate and unsorted read pairs

for inputbam in `ls ./*.bam`; do 
  prefixinput=`echo $inputbam | sed 's/.bam//i' | sed 's/\.\///i' `
  goodreads=`echo $prefixinput | sed 's/$/_goodreads/i' `
  samtools fastq $inputbam --threads 8 -N 2>error_samtools.tmp | skewer - --min 48 -n --quiet -1 --threads 8 --mode pe -x /opt/skewer-master/adapters/adapters.fa 2>error_skewer.tmp | perl /opt/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -min_len 48 -out_format 3 -out_good stdout -out_bad null -no_qual_header -min_qual_mean 10 -ns_max_p 5 -lc_method dust -lc_threshold 10 -trim_tail_left 5 -trim_tail_right 5 -trim_ns_left 1 -trim_ns_right 1 -trim_qual_left 10 -trim_qual_right 10 2> $goodreads.prinseq.log | gzip > $goodreads.fq.gz
done

# Fix pairs and save them in a R1, R2 and singletons separate files
for gzipped in `ls *_goodreads.fq.gz`; do 
  /scratch/lola/exception/single_paired_filter_fj.py $gzipped
done

