#!/bin/bash
### input is BAM with prealigned PAIRED-END reads to host genome file (human, mouse). For Single end reads should be different (ask Lola)

cd .../exception/HPGQ     # in this directory bam files with prealigned files are located. MAY THIS BE A VBLE TO SUBSTITUTE ALSO IN LINE 18, 25, 

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
  .../exception/single_paired_filter_fj.py $gzipped
done

# in the cluster (pero esto da igual... se puede hacer todo en la misma máq)
#!/bin/bash
### input is FASTQ ### Use directly bwa ALN without the bwa.kit. FASTQ is PE
# Usage: ./02-BWA-MEM-process_fastq_PE.sh fastq.fq.gz
cd /home/lalonso/projects/goldStd/goodreads/LPBQ
inputfq=$1
echo "........... processing file $inputfq .........."
prefixinput=`echo $inputfq | sed 's/_goodreads_r1.fq.gz//i'`
inputfq2=`echo $inputfq | sed 's/_r1.fq.gz/_r2.fq.gz/i'`
inputreads=`expr $(wc -l $inputfq | cut -f1 -d ' ') / 2`  # I will count 2, only r1 y r2
# Representative refSeq
bwa aln -t8 /home/lalonso/dbs/reprRefSeq/reprBact $inputfq > $prefixinput.r1.sai
bwa aln -t8 /home/lalonso/dbs/reprRefSeq/reprBact $inputfq2 > $prefixinput.r2.sai
bwa sampe /home/lalonso/dbs/reprRefSeq/reprBact $prefixinput.r1.sai $prefixinput.r2.sai $inputfq $inputfq2 | samtools view -1 - -h -f3 -F 1796 -q0 -b > $prefixinput.mapped.aln.bam  # Remove unmapped, not primary and QC fails (probably not working since we are in fastq)
samtools sort -o $prefixinput.sorted.bam $prefixinput.mapped.aln.bam
samtools flagstat $prefixinput.sorted.bam > $prefixinput.stats 
echo "original input fastq reads: $inputreads" >> $prefixinput.stats
# echo "non-aligned-to-human reads: $nonhumanreads" >> $prefixinput.stats
mv $prefixinput.stats /home/lalonso/projects/goldStd/results/reprBact
mv $prefixinput.sorted.bam /home/lalonso/projects/goldStd/results/reprBact/$prefixinput.bam
rm $prefixinput.mapped.aln.bam
rm $prefixinput.r1.sai
rm $prefixinput.r2.sai

# # # treatment of alignment bamfiles. This is common for databases, but has to be done separately ###
# Representative refSeq
cd /home/lalonso/projects/goldStd/results/reprBact/
outbamfile=$prefixinput.bam
samtools view $outbamfile | cut -f3 > "$prefixinput"all_appearances.txt
sort -u "$prefixinput"all_appearances.txt | grep -v -w '*' > "$prefixinput"unique_spp.txt
for sp in `cat "$prefixinput"unique_spp.txt `; do printf $sp"\t"; grep -w $sp -c "$prefixinput"all_appearances.txt; done > "$prefixinput"counts.txt
for sp in `cat "$prefixinput"unique_spp.txt `; do grep -w -m1 $sp /home/lalonso/dbs/reprRefSeq/reprBact_headers.tsv; done | cut -f2 > "$prefixinput"taxa.txt
paste "$prefixinput"counts.txt "$prefixinput"taxa.txt | sort -k2 -nr > "$prefixinput"_specific_sp.tsv
rm "$prefixinput"*.txt
samtools index $prefixinput.bam
samtools idxstats $prefixinput.bam | cut -f1,2,3  > "$prefixinput"_all_spp.tsv

#!/bin/bash
### we have a "$prefixinput"_all_spp.tsv output counts from samtools idxstats and have to merge them in a unique global table### 
# Usage:/home/lalonso/projects/goldStd/cleanReads/scripts/createGlobalTable.sh
# reprBact
cd /home/lalonso/projects/goldStd/results/reprBact/HPGQ  # LPBQ HPBQ HPGQ LSHPGQ
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
paste reprBact_counts.tmp /home/lalonso/dbs/reprRefSeq/reprBact_headers.tsv > reprBact_counts.tsv
header=`head -1 reprBact_counts.tsv`; header=`printf "$header%s\t%sTotal"`
awk 'BEGIN{FS="\t"; OFS="\t"} NR>1 {for(sample=3;sample<=NF-2;sample++) total+=$sample; print $0"\t"total; total=0}' reprBact_counts.tsv > reprBact_sorted.tmp
sort --numeric-sort --reverse --field-separator=$'\t' --key=$numbersamples reprBact_sorted.tmp > reprBact_counts_sorted.tsv
sed -i "1i $header" reprBact_counts_sorted.tsv
sed -i "$ d" reprBact_counts_sorted.tsv # remove last line that is always an asterisk
rm *.tmp


# BEDtools calculate coverage