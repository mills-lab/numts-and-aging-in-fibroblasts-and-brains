#October 2021
#This script describes steps to use 'Roadmap epigenome states' dataset to find overlaps with Numt positions in the genome. The end goal is to predict chromatin states that Numts may overlap with in human genome.
#Step1: Choose a bed file that is relevant to the sample/tissue-type for matching. 
#Step2: Import it with 'rtracklayer' as GRanges. Extract ranges from your vcf file and use "countOverlaps".
#Example dataset is bedfile 'E73' a 15-state epigenome data for DLPFC region. https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15

rm(list = ls())
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(MTseeker)
library(MTseekerData)
library(metablastr)
library(biomaRt)

#Define file paths for bed file as well as Numt vcf file.
dir<-"/filepath"
setwd("/filepath")

#Import bed file
e73<-import.bed("E073_15_coreMarks_dense.bed.gz", genome = "hg19")

#Import Numt call vcf file
vcf<- readVcf("filename.vcf")

#Filter calls by 'PASS' in the vcf file
vcf<- vcf[rowRanges(vcf)$FILTER == "PASS"]

#Define row ranges and sort, check for exact match
numt<-rowRanges(vcf)
numt<-sort(vcf) #vcf file
stopifnot(width(vcf) == 1)

#seqlevels should be matched between Numt and Bed file before counting overlaps
seqlevelsStyle(e73) <- "UCSC"
seqlevelsStyle(numt) <- "UCSC"
seqlevels(numt)<-paste0("chr", seqlevels(numt))
table(seqnames(e73))
table(seqnames(numt))

e73<- sort(e73)
e73l <- split(e73, e73$name)
e373l = as(e73l, "data.frame")
#optional to export data as csv
write.csv(e73l, "filename.csv")

#Function to count overlaps
c <- sapply(e73l, function (state) {
  return(sum(countOverlaps(numt, state)))
})
c
write.csv(c, "filename.csv")

#To perform binomial probability distribution to know enrichment of numts in chromatin states, test expected and observed state occurrences and their significance.
#if there was exactly 1 match per interval, then sum of all states should be equal to total numts in vcf file.

sum(c) == length(numt)
stateWidths <- sapply(e73l, function (state) {sum(width(state))})
stateProp <- stateWidths/sum(stateWidths)

#To calculate individual probabilities: pbinom(q=c["5_TxWk"] - 1, size=sum(c), prob=stateProp["5_TxWk"], lower.tail=FALSE)

for (i in 1:15) {
  pbinoml<-pbinom(q=c[i] - 1, size=sum(c), prob=stateProp[i], lower.tail=FALSE)
  write.table(pbinoml, "filename.csv", append=TRUE,sep=",",row.names=TRUE)
}

#*********END***************
