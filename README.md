# DRAC
Discarded Reads Alignment and Coverage

DRAC is conceived as a post-processing tool to be run on a NGS file (either fastq or SAM/BAM format) of previously discarded reads that could not be aligned to the host genome (it may be human or any other host). It is mostly devised to be run on a single organism tissue RNA-Seq (also Whole Exome Sequencing or Whole Genome Sequencing are supported) rather than a high biomass metagenomic design (i.e., it is not recommended for high biomass studies, e.g. feces or saliva). This is because the program that DRAC uses for identifying the microbes is the well know Burrows-Wheeler aligner (BWA), that it is not designed for that sort of job, besides of being slow due to the alignment process. Instead, its suitable to align a subset of Illumina paired-end (or single) sequencing data from previously aligned reads to the human (or host) reference genome, to complete bacterial genomes in the RefSeq database (as was performed in a seminal paper that identified bacteria that may be present in human tumor samples (DOI: 10.1186/s40168-016-0224-8 )

With BWA identifications, one cannot conclude whether the reads are true or false positives, DRAC pipeline afterwards will apply a method to ensure that counts belonging to a particular genome are well spread and covering a directly proportional part of it. 

DRAC is is a binner tool that takes a BAM file as input (with a minor change it can start from a non-host reads fastq file) that keeps both the host aligned reads and the unmapped (Discarded by the aligner, for whatever reason), that DRAC will classify by the bitwise flag (samtools). Afterwards, DRAC trims the putative remaining adapters, trims the low complexity sequences with the DUST algorithm, trims the low quality tails and in dismisses the short reads after this process, saving them in a fastq file. 

Since the remaining reads will be unsorted and one of the pairs members is sometimes lost, the pipeline fixes the pairs to keep only the paired reads. It will also save the singletons in case the user needs to analyze them in a further step.

BWA will align the unmapped, but sorted paired-end reads, against a reference database. In the article that we present here, the database is previously built (BWA index) with the bacterial genomes downloaded from RefSeq, and in particular, the set of Representative Genomes (one assembly per defined species selected as representative, with criteria explained in https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#representative_genomes) 

The alignment counts will be merged in a sample table (individual sample, or for a group of samples).

Afterwards, DRAC will use bedtools coverageBed for computing the coverage along the genomes represented in the sequence alignments of the resulting BAM file. The computed coverages will be merged in a sample table (individual sample, or for a group of samples).

The unique aspect of DRAC pipeline is that it will take into account the coverage of a particular bacterial genome and will evaluates the expected vs. the observed coverage.  DRAC assumes a constant (it may be the average, also calculated by the method) read length, and it also assumes that reads should be widespread or scattered along the genome (KrakenUniq follows a similar idea, https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0).

Following these principles, we developed when we know how long the reads are, and how much they should cover, we can calculate an internal score to discard those reads with low coverage (according to what it would be expected). By doing so, DRAC penalizes bacterial assignments with very poor coverage and will dismiss them, while it will retain the most reliable (scattered through the genome) ones, most probably, true positives.

The pipeline also performs some calculations, proportions, creates a global table in case of studying several samples, and in the last steps it will perform the counts normalization and statistical analysis, but these aspects are beyond the scope of the present article.

DRAC is not a metagenomics classifier, but instead a pipeline that was developed to scan Next-generation sequencing (NGS) files for putatively present microbiota or contaminant microbes.
