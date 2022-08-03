#!/bin/bash
### input is BAM 
### Use directly BWA on the BAM file
cd /home/lola/projects/goldStd/exception/LSHPGQ/
inputbam=$1
echo "........... processing file $inputbam .........."
prefixinput=`echo $inputbam | sed 's/.bam//i'`
inputreads=`samtools flagstat $inputbam | head -1 | cut -f1 -d ' '`
nonhumanbam="$prefixinput.non-human.bam"
samtools view -f 12 -F 1536 -b $inputbam > $nonhumanbam   # Keep unaligned (read unmapped + mate unmapped), remove QC fails or PCR dups
nonhumanreads=`samtools flagstat $nonhumanbam | head -1 | cut -f1 -d ' '`

# Alignment to NCBI RefSeq 
# Representative refSeq (1 microbe representative for its species) https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#representative_genomes
bwa aln -t4 /home/lola/dbs/reprRefSeq/reprBact -b1 $nonhumanbam > $prefixinput.R1.sai
bwa aln -t4 /home/lola/dbs/reprRefSeq/reprBact -b2 $nonhumanbam > $prefixinput.R2.sai
bwa sampe /home/lola/dbs/reprRefSeq/reprBact $prefixinput.R1.sai $prefixinput.R2.sai $nonhumanbam $nonhumanbam | samtools view -1 - -h -f 3 -F 1548 -q0 -b > $prefixinput.mapped.aln.bam # Keep paired and proper, remove unmapped & not primary (bitwise flag 1804). When we want to retain secondary alignments set the bitwise flag 1548
samtools sort -o $prefixinput.sorted.bam $prefixinput.mapped.aln.bam
samtools flagstat $prefixinput.sorted.bam > $prefixinput.stats
echo "original inputbam reads: $inputreads" >> $prefixinput.stats
echo "non-aligned-to-human reads: $nonhumanreads" >> $prefixinput.stats
mv $prefixinput.stats /home/lola/projects/goldStd/results/reprBact/
mv $prefixinput.sorted.bam /home/lola/projects/goldStd/results/reprBact/$prefixinput.bam
rm $prefixinput.mapped.aln.bam
rm $prefixinput.R1.sai
rm $prefixinput.R2.sai

# Extract stat info of alignment bamfiles. ###
# Representative refSeq
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

# convert to fasta to check homology in blast
# cd /home/lalonso/projects/goldStd/results/reprBact/
# samtools bam2fq $prefixinput.bam 2>bam2fq.stderr | awk 'NR%4==1 || NR%4==2' | sed 's/^@/>/g' > ./fasta2blast/$prefixinput.fasta
# rm bam2fq.stderr
