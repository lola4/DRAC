#!/bin/bash
# calculate the length of the reads that aligned to bacteria, because in some datasets they are trimmed very differently
# assuming HPGQ scenario described in the paper
PATHTOBAMFILES=/home/lola/projects/goldStd/results/reprBact/HPGQ
cd $PATHTOBAMFILES

rm lengthOfMappedReads.tsv  # in case there exists a previously generated file (previous tests)
touch lengthOfMappedReads.tsv
printf "sample\tlength\n" > lengthOfMappedReads.tsv
for bamfile in `ls *bam`; do 
  sample=`echo $bamfile | sed 's/.bam//i'`
  statsfile=`echo $bamfile | sed 's/bam/stats/i'`
  sumOfLengths=`samtools view $bamfile | cut -f10 | wc -c`  # counts one more character/seq (the tab?) so I have to subtract the total of summed seqs 
  mappedreads=`cat $statsfile | head -1 | cut -f1 -d ' '`
  average=`expr $(($sumOfLengths - $mappedreads)) / $mappedreads`
  printf $sample"\t"$average"\n" >> lengthOfMappedReads.tsv
done

# the output will be the nucleotide length (e.g. 84 bp/read, 125 bp/read) per each sample/library (usually it is the same and and integer number)


