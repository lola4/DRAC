#!/usr/bin/R
# From tables *_sorted.tsv (counts table and coverage table from the same dataset) this script will apply the penalization for genomes that are poorly covered

########################################
##### Read counts and coverage files ###
########################################
#  The unique aspect of DRAC pipeline is that it will take into account the coverage of a particular bacterial genome and it will evaluate the expected vs. the observed coverage.  
baseDir = "~/projects/LSHPGQ/"
setwd(paste0(baseDir,"analysis/"))
database = "reprBact"  

# DRAC assumes a constant read length (the average length, also calculated by the method, in the lengthOfMappedReads.tsv file), and it also assumes that reads should be widespread or scattered along the genome for them to be considered reliable (KrakenUniq follows a similar idea, and so does SLIMM, discarding unlikely genomes based on coverage landscape). This way, those methods eliminate genomes that have a stack of reads only in few places across their genomes, which could be a result of sequencing artifact, spike-ins, or a conserved region among distant relatives.
lengthOfMappedReads = 48 # may be any number (100bp, 150bp, also non-integers)

counts_table <- as.data.frame(read.delim(file=paste0(baseDir,database,"/LSHPGQ/",database,"_counts_sorted.tsv"), header=FALSE, sep="\t", comment.char="", stringsAsFactors=FALSE))
rownames(counts_table) = counts_table$V1

# calculate abundances (reads aligned to bacterial species normalized to alignedReads)
abundances = subset(counts_table[-c(1),], select=-c(V1,V2,V24,V25,V26))
abundances <- as.data.frame(lapply(abundances, as.integer))
rownames(abundances) = rownames(counts_table[-c(1),])
abundances = 1*t(apply(abundances, 1, `/`,1)) # Transform into dataframe instead of list. We have to keep raw counts for coverage

### ### #################################################### ###
### ### cov. penalization (avoid piled up reads) ### ### ### ###
### ### ### ################################################ ###
coverage_table <- as.data.frame(read.delim(file=paste0(baseDir,database,"/LSHPGQ/",database,"_coverage_sorted.tsv"), header=FALSE, sep="\t", comment.char="", stringsAsFactors=FALSE))
rownames(coverage_table) = coverage_table$V1
covtable = subset(coverage_table[-c(1),], select=-c(V1,V2,V24,V25,V26)) # This is to remove the text cells: taxon name, etc
identical(colnames(covtable), colnames(abundances))  # check that coverage table is in the same order of counts_table (abundances)

coverAsCounts <- covtable[rownames(abundances),]
covtable <- as.data.frame(lapply(coverAsCounts, as.integer))
covtable = 1*t(apply(covtable, 1, `*`,1)) # do not normalize but transform into dataframe instead of list
rownames(covtable) = rownames(coverAsCounts)
identical(rownames(covtable), rownames(abundances))

checkCovCounts <- function(subject) { # also check that the abundances matrix has value nonzero in the same places than the coverage
  identical(covtable[,subject] == 0 , abundances[,subject] == 0 )
}
table(t(sapply(colnames(covtable), checkCovCounts)))  # sapply retrieves a transposed matrix
rm(list=c("coverAsCounts"))

### ### ### ### ### ### #################################################### ###
### ### ### Extract TAXA with POSITIVE counts only & NORMALIZE COVERAGE  ### ###
### ### ### ### ### ### #################################################### ###
# exclude absent species from the dataset and create a dataframe with positive_counts species only
positive_counts <- abundances # 4331
positive_counts <- positive_counts[rowSums(positive_counts) > 0,]
identical(colnames(abundances), colnames(positive_counts))
# do the same for the coverages (should be the same number of positive taxa)
positive_coverages <- covtable # 4331
positive_coverages <- positive_coverages[rowSums(positive_coverages) > 0,]
identical(rownames(positive_coverages), rownames(positive_counts))

# calculate the (max) observed_coverages_vs_theoretical ratio. 
obsExpWeights = matrix(nrow=length(rownames(positive_counts)), ncol=dim(positive_counts)[2])
colnames(obsExpWeights) = colnames(positive_counts)
rownames(obsExpWeights) = rownames(positive_counts)
matObsExpCov <- function(taxon){ 
  obsExpWeights[taxon,] = positive_coverages[taxon,]/(lengthOfMappedReads*positive_counts[taxon,])  # 84bp, 50bp
}
obsExpMatrix = t(sapply(rownames(positive_counts), matObsExpCov))
obsExpMatrix[is.nan(obsExpMatrix)] <- 0
dim(obsExpMatrix)


### ### ### ### ### ### #################################################### ###
### ###   Filter spurious TAXA with matObsExpCov 0-1 & NORMALIZE COUNTS  ### ###
### ### ### ### ### ### #################################################### ###
# We apply the DRAC algorithm: when we know how long the reads are, and how much they should cover if they are fairly spread (read length x number of reads), we can calculate an internal score to discard those reads with low coverage (according to what it would be expected). By doing so, DRAC penalizes bacterial assignments with very poor coverage and dismisses them, while it will retain the most reliable (scattered through the genome) ones, that most probably are true positives.
# we keep backup of several files in case we want to go back to the original matrices
positive_counts_original = positive_counts
positive_coverages_original = positive_coverages
obsExpMatrix_original = obsExpMatrix

# filter taxa that are under 50% of the expected maximum theoretical coverage. Flatten obsExpMatrix: less than 0.5 will become 0, more than 0.5 of exp_cov become 1. To let counts be counts and no other strange measure unit
thresholdCov = 0.1
obsExpMatrix[obsExpMatrix < thresholdCov] = 0
obsExpMatrix[obsExpMatrix >= thresholdCov] = 1
table(obsExpMatrix)
# multiply the matrix
identical(rownames(positive_counts),rownames(obsExpMatrix))
positive_counts = obsExpMatrix * positive_counts

positive_counts = 1*t(apply(positive_counts, 1, `/`,1))  # use Raw counts (collapse below).

# Remove bacteria that have become 0 in all samples to avoid using them in upstream analysis
positive_counts <- positive_counts[rowSums(positive_counts) > 0,]  # from 4331 initial taxa to a much lower number

write.table(positive_counts, file="reliable_counts.tsv", sep="\t", header=TRUE)

# In the last stages DRAC will perform the counts normalization and statistical analysis, but these aspects are beyond the scope of the present article.

