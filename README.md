# Associations of nuclear mitochondrial DNA insertions and human lifespan in aging fibroblasts and in the human brain

## Variant calling pipelines

Requires: 
```
exonerate/2.4.0 (for Dinumt) 
vcftools/0.1.15 (for Dinumt) 
bcftools/1.5
Dinumt/0.0.23
Delly2/0.8.5
bowtie2/2.1.0 (for MELT) 
MELT/2.1.4
```

For non-reference Numt calling:
```
${DINUMT}/dinumt.pl \
--mask_filename=${DINUMT}/dinumt/refNumts.bed \
--input_filename=${DATA}/bam/hfb12-1_S1.bam \
--reference=.${REF}/hg19.fa \
--min_reads_cluster=1 \
--include_mask \
--mt_names=chrM \
--output_filename=hfb12-1_S1.bam.final.vcf \
--prefix=hfb12-1_S1.bam.final \
--len_cluster_include=564 \
--len_cluster_link=1128 \
--insert_size=306 \
--max_read_cov=100 \
--output_support \
--support_filename=hfb12-1_S1.bam.final_support.sam \
--ucsc
```

For merging non-reference Numt VCFs into one file:
```
grep ^# hfb12-1_S1.bam.final.vcf  > header.txt 
cat *vcf | grep -v ^# | vcf-sort -c | perl ${DINUMT}/clusterNumtsVcf.pl --samtools=samtools --reference=${REF} > data.txt 
cat header.txt data.txt > merged.vcf
```

For non-reference SV calling:
```
${DELLY}/src/delly call -g ${REF}/hg19.fa -x ${DELLY}/excludeTemplates/human.hg19.excl.tsv ${DATA}/hfb12-1_S1.bam -o hfb12-1_S1.bam.delly.bcf
bcftools view hfb12-1_S1.bam.delly.bcf > hfb12-1_S1.bam.delly.vcf
```

For non-reference MEI calling:
```
java -jar ../MELTv2.1.4/MELT.jar Single -bamfile ${DATA}/hFB7-8_S1.bam -w ../hFB7-8_S1.bam.melt.0205 -t ../MELTv2.1.4/me_refs/1KGP_Hg19/LINE1_MELT.zip -h ${REF}/hg19.fa -n ../MELTv2.1.4/add_bed_files/1KGP_Hg19/hg19.genes.bed -c 21
java -jar ../MELTv2.1.4/MELT.jar Single -bamfile ${DATA}/hFB7-8_S1.bam -w ../hFB7-8_S1.bam.melt.0205 -t ../MELTv2.1.4/me_refs/1KGP_Hg19/ALU_MELT.zip -h ${REF}/hg19.fa -n ../MELTv2.1.4/add_bed_files/1KGP_Hg19/hg19.genes.bed -c 21
java -jar ../MELTv2.1.4/MELT.jar Single -bamfile ${DATA}/hFB7-8_S1.bam -w ../hFB7-8_S1.bam.melt.0205 -t ../MELTv2.1.4/me_refs/1KGP_Hg19/SVA_MELT.zip -h ${REF}/hg19.fa -n ../MELTv2.1.4/add_bed_files/1KGP_Hg19/hg19.genes.bed -c 21
```

## Epigenetic analysis

This script describes steps to use 'Roadmap epigenome states' dataset to find overlaps with Numt positions in the genome. The end goal is to predict chromatin states that Numts may overlap with in human genome.

Step1: Choose a bed file that is relevant to the sample/tissue-type for matching. 

Step2: Import it with 'rtracklayer' as GRanges. Extract ranges from your vcf file and use "countOverlaps".

Example dataset is bedfile 'E73' a 15-state epigenome data for DLPFC region. https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15

Requires: 
```
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(MTseeker)
library(MTseekerData)
library(metablastr)
library(biomaRt)
```

Run
```
/lib/R_script_epigenome.R
```

## 

Requires: 
```

```

Run
```

```


## Citation


For Dinumt:
* Dayama, Gargi, Sarah B Emery, Jeffrey M Kidd, and Ryan E. Mills. 2014. [The genomic landscape of polymorphic human nuclear mitochondrial insertions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4227756/pdf/gku1038.pdf),
Nucleic Acids Research, 2014, gku1038, `https://doi.org/10.1093/nar/gku1038`


## Contact:

arthurz@umich.edu or https://github.com/WeichenZhou

kalpita.karan@gmail.com or https://github.com/kalpita23