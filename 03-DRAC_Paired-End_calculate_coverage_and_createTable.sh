#!/bin/bash
### # BEDtools calculate coverage
# no input needed, except for the directory where the BAM files are located $PATHTOBAMFILES. For example # cd $PATHTOBAMFILES/coverage
# assuming HPGQ scenario described in the paper
PATHTOBAMFILES=/home/lola/projects/goldStd/results/reprBact/HPGQ
cd $PATHTOBAMFILES.
mkdir coverage
cd /home/lola/projects/goldStd/results/reprBact/HPGQ/coverage
mkdir genomeCovBed

# calculate genomeCoverageBed for every bamfile [13 min]. This command is the same as bedtools coverage. https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
for bamfile in `ls $PATHTOBAMFILES/*.bam`; do covFile=`echo $bamfile | sed 's/^.*HPGQ\///' | sed 's/.bam$/.cov/' `; genomeCoverageBed -ibam $bamfile -bga -split > $PATHTOBAMFILES/coverage/$covFile; done
mv *.cov ./coverage/genomeCovBed/

# sum the covered bases along the genome (these covered bases will be used to penalize the reads that pile up in one single location) 
for covfile in `ls $PATHTOBAMFILES/coverage/genomeCovBed/*.cov`; do covbpFile=`echo $covfile | sed 's/^.*genomeCovBed\///' `; awk 'BEGIN{FS="\t"; OFS="\t"} {if ($4>0) {$4=1}; print $0}' $covfile | awk 'BEGIN{FS="\t"; OFS="\t"} {result=($3-$2)*$4} {print $1,result}' | awk 'BEGIN{FS="\t"; OFS="\t"} {taxa[$1]+=$2} END{for(taxon in taxa) print taxon,taxa[taxon] } ' > $covbpFile; done

# now we have to sum up all of this in one table
templatefile=`ls ..../results/reprBact/HPGQ/*_all_spp.tsv | head -1`
cut -f1,2 $templatefile > reprBact_counts.tsv
# we have to remove the last line that is always an asterisk (belonging to the unmapped reads)
mv reprBact_counts.tsv reprBact_countswithtrail.tsv
head -n-1  reprBact_countswithtrail.tsv > reprBact_counts.tsv

# create accumulated presences/absences file
rm accum_presences.tsv
touch accum_presences.tmp
touch accum_presences.tsv
for taxon in `cut -f1 reprBact_counts.tsv`; do
  printf $taxon"\t" > presence.tmp
  grep -w $taxon *.cov -c -h | tr '\n' '\t' >> presence.tmp
  echo  >> presence.tmp
  cat accum_presences.tsv presence.tmp > accum_presences.tmp
  cp accum_presences.tmp accum_presences.tsv
done
rm *.tmp

# create header for global coverage file
printf 'RefSeq_ID\tSeq_Length\t' > header.txt
numbersamples=`expr $(ls *.cov | wc -l) + 5`
for covFile in `ls $PATHTOBAMFILES/coverage/*.cov`; do
  sample=`echo $covFile | sed 's/$PATHTOBAMFILES\/coverage\///i' | sed 's/\.cov//i'`;
  printf $sample"\t" >> header.txt
done
header=`head -1 header.txt`; header=`printf "$header%sRefSeq_ID\t%sdefinition_NCBI_RefSeq\t%sTotal"`

# create global coverage file
rm coverages_reprBact.tsv
touch coverages_reprBact.tsv
for taxon in `cut -f1 reprBact_counts.tsv `; do
  printf $taxon"\t" >> coverages_reprBact.tsv
  for covFile in `ls $PATHTOBAMFILES/coverage/*.cov`; do
    outputGrep=`grep -w -P $taxon"\t" -m1 $covFile | cut -f2`
    if [ -z "$outputGrep" ]; then
       cellContent=0
    else 
       cellContent=$outputGrep
    fi 
    printf $cellContent"\t" >> coverages_reprBact.tsv
  done
  printf "\n" >> coverages_reprBact.tsv
done

# create a nice table also sorted by top most covered species genomes 
paste reprBact_counts.tsv coverages_reprBact.tsv > coverages_reprBact0.tmp
awk 'BEGIN {FS="\t"} {if ($1 != $3) print $1,$3}' coverages_reprBact0.tmp  # check that bacteria are equal after the paste
cut -f3 --complement coverages_reprBact0.tmp > coverages_reprBact1.tmp
for sp in `cut -f1 coverages_reprBact1.tmp`; do grep -w -m1 $sp /local/microbiome/dbs/reprBact/reprBact_headers.tsv; done > bacteriataxa.tmp
paste coverages_reprBact1.tmp bacteriataxa.tmp > coverages_reprBact2.tmp
awk 'BEGIN {FS="\t"} {if ($1 != $25) print $1,$25}' coverages_reprBact2.tmp  # check that bacteria are equal after the paste
cut -f24 --complement coverages_reprBact2.tmp > coverages_reprBact3.tmp
awk 'BEGIN{FS="\t"; OFS="\t"} {for(sample=3;sample<=NF-2;sample++) total+=$sample; print $0"\t"total; total=0}' coverages_reprBact3.tmp > coverages_reprBact4.tmp
sort --numeric-sort --reverse --field-separator=$'\t' --key=$numbersamples coverages_reprBact4.tmp > reprBact_coverage_sorted.tsv
sed -i "1i $header" reprBact_coverage_sorted.tsv

rm *tmp
rm header.txt

mv reprBact_coverage_sorted.tsv $PATHTOBAMFILES/

# the table reprBact_coverage_sorted.tsv now contains all the reads covered per sample, for each genome/contig, with species in rows and samples in columns. It has the same format than the reprBact_counts_sorted.tsv table with the numer of read counts aligned, with the samples and species in the exact same order.
# The column identifiers will be the names of the input BAM files, and the row identifiers will be the NCBI accession numbers
# the table reprBact_coverage_sorted.tsv, along with the reprBact_counts_sorted.tsv, will be input for the last script of the DRAC pipeline:  05-apply_Cov_penalization.R
