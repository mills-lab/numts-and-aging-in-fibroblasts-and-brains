# Associations of nuclear mitochondrial DNA insertions and human lifespan in aging fibroblasts and in the human brain

## 

Requires: 
```

```

Run
```

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