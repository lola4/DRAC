#!/bin/bash
# we have a "$prefixinput"_all_spp.tsv output counts from samtools idxstats and have to merge them in a unique global table that will summarize the reads aligned to the set of Representative Genomes (reprBact) for all samples.

cd /home/lola/projects/goldStd/results/reprBact/HPGQ  # LPBQ HPBQ HPGQ LSHPGQ # different scenarios from paper
templatefile=`ls *_all_spp.tsv | head -1`
cut -f1,2 $templatefile > reprBact_counts.tsv
sed -i '1i RefSeq_ID\tSeq_Length' reprBact_counts.tsv
numbersamples=`expr $(ls *_all_spp.tsv | wc -l) + 5`
for countsfile in `ls *_all_spp.tsv`; do
  sample=`echo $countsfile | sed 's/_all_spp.tsv//i'`;
  echo $sample > $sample.tmp
  cut -f3 $countsfile >> $sample.tmp
  paste reprBact_counts.tsv $sample.tmp > reprBact_counts.tmp
  cp reprBact_counts.tmp reprBact_counts.tsv
done
paste reprBact_counts.tmp /home/lola/dbs/reprRefSeq/reprBact_headers.tsv > reprBact_counts.tsv
header=`head -1 reprBact_counts.tsv`; header=`printf "$header%s\t%sTotal"`
awk 'BEGIN{FS="\t"; OFS="\t"} NR>1 {for(sample=3;sample<=NF-2;sample++) total+=$sample; print $0"\t"total; total=0}' reprBact_counts.tsv > reprBact_sorted.tmp
sort --numeric-sort --reverse --field-separator=$'\t' --key=$numbersamples reprBact_sorted.tmp > reprBact_counts_sorted.tsv
sed -i "1i $header" reprBact_counts_sorted.tsv
sed -i "$ d" reprBact_counts_sorted.tsv # remove last line that is always an asterisk
rm *.tmp

# the table reprBact_counts_sorted.tsv now contains all the reads aligned per sample, with species in rows and samples in columns. The column identifiers will be the names of the input BAM files, and the row identifiers will be the NCBI accession numbers
# the table reprBact_counts_sorted.tsv will be input for the last script of the DRAC pipeline:  05-apply_Cov_penalization.R
