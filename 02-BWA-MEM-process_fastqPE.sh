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

# # convert to fasta to check homology in blast
# cd /home/lalonso/projects/goldStd/results/reprBact/
# samtools bam2fq $prefixinput.bam 2>bam2fq.stderr | awk 'NR%4==1 || NR%4==2' | sed 's/^@/>/g' > ./fasta2blast/$prefixinput.fasta
# rm bam2fq.stderr
