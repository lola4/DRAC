#!/bin/bash
### Remove putative remaining adapters, remove Low Complexity reads, trim Low Quality ends, remove too short reads
# input is a BAM file with the non-human (non-host) reads that have been discarded by the aligner (SAM flags 4, or combinations, like SAM flag 13).
# https://broadinstitute.github.io/picard/explain-flags.html

# Ths script performs the following actions on a BAM files (already mapped to a host genome) (can be applied recursively on a set of BAM files in a directory)
# 1- convert bam to fastq
# 2- trim adapters or remove reads with adapters with skewer.
# 3- remove low complexity and mean low qual with prinseq
# 4- save good reads in fastq.gz that will contain separate and unsorted read pairs
# 5- fix pairs
# 6- aligns the reads to the set of Representative Genomes (downloaded previously from https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#representative_genomes)

# Sequence cleaning
inputbam=$1
prefixinput=`echo $inputbam | sed 's/.bam//i' | sed 's/\.\///i' `
goodreads=`echo $prefixinput | sed 's/$/_goodreads/i' `
samtools fastq $inputbam --threads 8 -N 2>error_samtools.tmp | skewer - --min 48 -n --quiet -1 --threads 8 --mode pe -x /opt/skewer-master/adapters/adapters.fa 2>error_skewer.tmp | perl /opt/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -min_len 48 -out_format 3 -out_good stdout -out_bad null -no_qual_header -min_qual_mean 10 -ns_max_p 5 -lc_method dust -lc_threshold 10 -trim_tail_left 5 -trim_tail_right 5 -trim_ns_left 1 -trim_ns_right 1 -trim_qual_left 10 -trim_qual_right 10 2> $goodreads.prinseq.log | gzip > $goodreads.fq.gz

# Fix pairs and save them in a R1, R2 and singletons separate files
gzipped=$goodreads.fq.gz
single_paired_filter_fj.py $gzipped

# Sequence alignment + statistics
inputfq=$gzipped
echo "........... processing file $inputfq .........."
prefixinput=`echo $inputfq | sed 's/_goodreads_r1.fq.gz//i'`
inputfq2=`echo $inputfq | sed 's/_r1.fq.gz/_r2.fq.gz/i'`
inputreads=`expr $(wc -l $inputfq | cut -f1 -d ' ') / 2`  
# Representative refSeq
bwa aln -t8 /home/lola/dbs/reprRefSeq/reprBact $inputfq > $prefixinput.r1.sai
bwa aln -t8 /home/lola/dbs/reprRefSeq/reprBact $inputfq2 > $prefixinput.r2.sai
bwa sampe /home/lola/dbs/reprRefSeq/reprBact $prefixinput.r1.sai $prefixinput.r2.sai $inputfq $inputfq2 | samtools view -1 - -h -f3 -F 1796 -q0 -b > $prefixinput.mapped.aln.bam  # Remove unmapped, not primary and QC fails
samtools sort -o $prefixinput.sorted.bam $prefixinput.mapped.aln.bam
samtools flagstat $prefixinput.sorted.bam > $prefixinput.stats
echo "non-human input fastq reads: $inputreads" >> $prefixinput.stats
mv $prefixinput.stats /home/lola/projects/goldStd/results/reprBact
mv $prefixinput.sorted.bam /home/lola/projects/goldStd/results/reprBact/$prefixinput.bam
rm $prefixinput.mapped.aln.bam
rm $prefixinput.r1.sai
rm $prefixinput.r2.sai

# Extract how many reads for each Representative refSeq genome
cd /home/lola/projects/goldStd/results/reprBact/
outbamfile=$prefixinput.bam
samtools view $outbamfile | cut -f3 > "$prefixinput"all_appearances.txt
sort -u "$prefixinput"all_appearances.txt | grep -v -w '*' > "$prefixinput"unique_spp.txt
for sp in `cat "$prefixinput"unique_spp.txt `; do printf $sp"\t"; grep -w $sp -c "$prefixinput"all_appearances.txt; done > "$prefixinput"counts.txt
for sp in `cat "$prefixinput"unique_spp.txt `; do grep -w -m1 $sp /home/lola/dbs/reprRefSeq/reprBact_headers.tsv; done | cut -f2 > "$prefixinput"taxa.txt
paste "$prefixinput"counts.txt "$prefixinput"taxa.txt | sort -k2 -nr > "$prefixinput"_specific_sp.tsv
rm "$prefixinput"*.txt
samtools index $prefixinput.bam
samtools idxstats $prefixinput.bam | cut -f1,2,3  > "$prefixinput"_all_spp.tsv

# run script 02-DRAC_Paired-End_create_table.sh
