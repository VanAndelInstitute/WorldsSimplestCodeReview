#Did not run these.
---
title: "Code review: TET2 and hypermethylation"
author: "Tim Triche"
date: "November 22nd, 2021"
output: 
  html_document:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{TET2}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


#Start code running here!!!

# First:Installation

Install the WorldsSimplestCodeReview package, if you haven't. 

#Install these packages and run code w/o hashtag
```r
#install.packages("remotes")
#install.packages("BiocManager")
#BiocManager::install("VanAndelInstitute/WorldsSimplestCodeReview")
library(knitr)
```

#This code reads an error, requires additional? When ran, a '+' sign is added to next line. Excluded this code as a result.To remove lus sign hit "ESC" 
To extract just the R code, you can use knitr::knit(input, tangle=TRUE):



# Introduction

Long before any of you were born, back in 2010, an exciting paper came out 
which purported to show that _IDH1_, _IDH2_, and _TET2_ mutations shared a 
phenotype of hypermethylation owing to loss of 5-hydroxymethylcytosine. The 
details can be found in [the paper](https://doi.org/10.1016/j.ccr.2010.11.015), 
which is indeed a landmark. Nevertheless, some fine details of the work seemed
to disagree with the results of other cohorts when replication was attempted.

![The money shot](figure/TET2.png)

Some of you who have seen volcano plots before can guess where this is going.

# The data


```r

library(limma)


library(GEOquery)

#DName not found, showed up as warning message
if (!exists("DNAme")) data(DNAme)


if (FALSE) { # this takes about 5 minutes:

  # needed to fetch data
  library(GEOquery) 
  install
  MSK_HOVON <- getGEO("GSE24505")

  # skip the expression data:
  platform <- sapply(MSK_HOVON, annotation)
  methylation <- which(platform == "GPL6604")
  DNAme <- MSK_HOVON[[methylation]] # GPL6604, HG17_HELP_PROMOTER 
  DNAme$male <-ifelse(DNAme$characteristics_ch1=="sex (male.1_female.2): 1",1,0)
  DNAme$TET2 <- ifelse(DNAme$characteristics_ch1.7 == "tet2: WT", 0, 1)
  DNAme$IDH <- ifelse(DNAme$characteristics_ch1.8 == "idh1.idh2: WT", 0, 1)
  DNAme$purity <- as.integer(DNAme$"bm_%blasts:ch1") / 100
  save(DNAme, file="../data/DNAme.rda")

}

# how many probes, how many patients?
dim(DNAme)
#> Features  Samples 
#>    25626      394
# Features  Samples
#    25626      394
```

### Some contrasts

Is it the case that TET2, IDH1, and IDH2 mutations are exclusive?


```r

# always plot your data
library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.8.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
medata <- pData(DNAme)[, c("TET2", "IDH")]
mutations <- t(as.matrix(medata))
Heatmap(mutations,
    col = c("lightgray", "darkred"), name = "mutant", column_km = 4,
    column_names_gp = gpar(fontsize = 7)
)
```

![](TET2_files/figure-html/heatmap-1.png)<!-- -->

```r

# get any samples that have both genes mutated
(non_exclusive  <- medata %>% dplyr::filter(TET2 == 1 & IDH == 1))
#>           TET2 IDH
#> GSM604380    1   1
medata %>% table()
#>     IDH
#> TET2   0   1
#>    0 302  50
#>    1  41   1
```
❗ _GSM604380 is not exclusive for the mutation._ One of the results from the
paper is that the mutations are mutually exclusive. The data from this code
does not match what is found in their results.

Do we see genome-wide hypermethylation from TET2 mutations? 


```r

# model TET2 and IDH1/2 mutant related hypermethylation
# note: there are plenty of confounders (pb%, bm%, wbc) that could be included
library(limma) 

# simplest design: model gene expression with IDH and TET2 mutation status and 
# pull top ranked significant genes from liner model fit, get probes for each
# gene
design1 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 ))
fit1 <- eBayes(lmFit(exprs(DNAme), design1))
(IDH_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("IDH", colnames(design1)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
#> [1] 6513
# 6513 probes for IDH

(TET_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("TET2", colnames(design1)),
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
#> [1] 6
# 6 probes for TET2

# control for sex
design2 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 + male ))
fit2 <- eBayes(lmFit(exprs(DNAme), design2))
(IDH_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                           coef=grep("IDH", colnames(design2)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
#> [1] 6651
# 6651 probes for IDH 

(TET2_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                            coef=grep("TET", colnames(design2)),
                                            p.value=0.05, # change if you like 
                                            number=Inf)))
#> [1] 7
# 7 probes for TET2

# control for blast count
design3 <- with(pData(DNAme), model.matrix( ~ IDH:purity + TET2:purity))
fit3 <- eBayes(lmFit(exprs(DNAme)[, as.integer(rownames(design3))], design3))

(IDH_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                           coef=grep("IDH", colnames(design3)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
#> [1] 7450
# 7450 probes for IDH:purity

(TET2_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                            coef=grep("TET", colnames(design3)),
                                            p.value=0.05, # change if you like 
                                            number=Inf)))
#> [1] 10
# 10 probes for TET2:purity
```
