---
title: "Code review: TET2 and hypermethylation"
author: "VAIGS Intro to Experimental Design, Class of 2026"
date: "November 22nd, 2021"
output: 
  html_document:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{TET2}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



# Installation

Install the WorldsSimplestCodeReview package, if you haven't. 
#This chunk was not necessary to run for me.-Riley


```r
#install.packages("remotes")
#install.packages("BiocManager")
#BiocManager::install("VanAndelInstitute/WorldsSimplestCodeReview")
library("WorldsSimplestCodeReview")
```

To extract just the R code, you can use `knitr::knit(input, tangle=TRUE)`.

# Introduction

Long before any of you were born, back in 2010, an exciting paper came out 
which purported to show that _IDH1_, _IDH2_, and _TET2_ mutations shared a 
phenotype of hypermethylation owing to loss of 5-hydroxymethylcytosine. The 
details can be found in [the paper](https://doi.org/10.1016/j.ccr.2010.11.015), 
which is indeed a landmark. Nevertheless, some fine details of the work seemed
to disagree with the results of other cohorts when replication was attempted.

![The money shot](figure/TET2.png) 
# This is Fig6 from the paper, suggesting TET2 mutant AML is associated with a hypermethylation p/type. Volcano plot (D) shows skew towards higher amount fo methylation in TET2 mutant compared to both TET2 and IDH WT. NB IDH and TET2 mutations mutually exclusive in AML.

Some of you who have seen volcano plots before can guess where this is going.

# The data
Loading the data from the paper. This requred loading packages limma and GEOquery

```r
#Load a package that can be used for Data analysis, linear models, and differential expression for microarray data. This is part of the BioConductor package.
library(limma)
#Load GEOquery package, that allows you to get data from the NCBI Gene Expression Omnibus
library(GEOquery)
#> Loading required package: Biobase
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following object is masked from 'package:WorldsSimplestCodeReview':
#> 
#>     plotMA
#> The following object is masked from 'package:limma':
#> 
#>     plotMA
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:dplyr':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:base':
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted,
#>     lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
#>     pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
#>     tapply, union, unique, unsplit, which.max, which.min
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#Load the paper's dataset specifically, if it isn't already; but just the parts you want
if (!exists("DNAme")) data(DNAme)

if (FALSE) { # this takes about 5 minutes:
#This didn't actually take 5 minutes for me so I hope it hasn't missed anything! Updated none
  # needed to fetch data
  library(GEOquery) 
  MSK_HOVON <- getGEO("GSE24505") #pulls GSE24505 dataset from GEOquery and labels it MSK_HOVON

  # skip the expression data:
  platform <- sapply(MSK_HOVON, annotation) #Using MSK_HOVON, pull datapoints with annotation into "platform"
  methylation <- which(platform == "GPL6604") #Select those datapoints within "platform" for which values are "GPL6604"- methylation data
  DNAme <- MSK_HOVON[[methylation]] # GPL6604, HG17_HELP_PROMOTER. Create new dataset "DNAme" with just methylation data
  DNAme$male <-ifelse(DNAme$characteristics_ch1=="sex (male.1_female.2): 1",1,0) #Variable w/in DNAme for sex
  DNAme$TET2 <- ifelse(DNAme$characteristics_ch1.7 == "tet2: WT", 0, 1) #Variable w/in DNAme for TET2 mutation/ WT
  DNAme$IDH <- ifelse(DNAme$characteristics_ch1.8 == "idh1.idh2: WT", 0, 1) #Variable w/in DNAme for IDH mutation/WT
  DNAme$purity <- as.integer(DNAme$"bm_%blasts:ch1") / 100 #Variable w/in DNAme for % purity 
  save(DNAme, file="../data/DNAme.rda") #Save DNAme dataset in local files

}


# how many probes, how many patients?
dim(DNAme) #dimensions of new DNAme dataset; features == probes, samples == patients
#> Features  Samples 
#>    25626      394
# Features  Samples
#    25626      394

#SO, we have 25626 probes and 394 patients/samples
```

### Some contrasts

Is it the case that TET2, IDH1, and IDH2 mutations are exclusive?


```r
# always plot your data

#this is the type of plot we are using...
library(ComplexHeatmap)

#make a matrix of methylation data, coded as a binary
mutations <- t(as.matrix(pData(DNAme)[, c("TET2", "IDH")]))
#plot the heatmap, I decreased font size to 4 so that single double mutant name was visible
  #I like the color purple better than red
Heatmap(mutations, col=c("lightgray","purple"), name="mutant", column_km=4,row_km = 1,
        column_names_gp = gpar(fontsize = 4))
```

![plot of chunk heatmap](figure/heatmap-1.png)

```r

#Group 1 has no mutations, 2 has only IDH muations, 3 has only TET2 mutations, 4 has both TET2 and IDH mutations
```
Group 4 could have been easily missed if the data was not plotted!
GSM604380 has both IDH and TET2 mutations.

Do we see genome-wide hypermethylation from TET2 mutations? 
  # No?!


```r
#This code is looking at methylation probes based on three different designs
  #Design 1 controls for no extra variables
  #Design 2 controls for sex
  #Design 3 controls for purity

# model TET2 and IDH1/2 mutant related hypermethylation
# note: there are plenty of confounders (pb%, bm%, wbc) that could be included
library(limma) 



# simplest design: model gene expression with IDH and TET2 mutation status and 
# pull top ranked significant genes from liner model fit, get probes for each
# gene

design1 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 ))
  #then use this matrix to run bayesian analysis
fit1 <- eBayes(lmFit(exprs(DNAme), design1))
  # & now get number of sig (p<0.05) probes, split by IDH vs TET2

(IDH_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("IDH", colnames(design1)), 
                                           p.value=0.05, 
                                           number=Inf)))
#> [1] 6513
# 6513 probes methylated at a significance value of 0.05 for IDH

(TET_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("TET2", colnames(design1)),
                                           p.value=0.05, 
                                           number=Inf)))
#> [1] 6
# 6 probes methylated at a significance value of 0.05 for TET2

# control for sex
design2 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 + male )) #table with binary coding as design 1, but with addition of sex variable where 1==male, 0==female
fit2 <- eBayes(lmFit(exprs(DNAme), design2)) #second Bayesian model with design 2 ie including sex confounder

(IDH_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                           coef=grep("IDH", colnames(design2)), 
                                           p.value=0.05, 
                                           number=Inf)))
#> [1] 6651
# 6651 probes methylated at a significance value of 0.05 for IDH 
  # Higher number of significant probes for IDH when including sex variable

(TET2_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                            coef=grep("TET", colnames(design2)),
                                            p.value=0.05, # change if you like 
                                            number=Inf)))
#> [1] 7
# 7 probes methylated at a significance value of 0.05 for TET2
  # Slightly higher (+1) significant number of probes for TET2 when including sex variable

# control for blast count
design3 <- with(pData(DNAme), model.matrix( ~ IDH:purity + TET2:purity)) #Table looking at %purity of sample for both IDH and TET2 on continuous scale from 0-1
fit3 <- eBayes(lmFit(exprs(DNAme)[, as.integer(rownames(design3))], design3)) # Bayesian analysis for purity 


(IDH_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                           coef=grep("IDH", colnames(design3)), 
                                           p.value=0.05, 
                                           number=Inf)))
#> [1] 7450
# 7450 probes methylated at a significance value of 0.05 for IDH:purity
  # Higher number of significant probes when fitting for % DNA purity of sample for IDH
(TET2_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                            coef=grep("TET", colnames(design3)),
                                            p.value=0.05, 
                                            number=Inf)))
#> [1] 10
# 10 probes methylated at a significance value of 0.05 for TET2:purity
  #Higher number of significant probes when fitting for % DNA purity of sample for TET2
```
# Conclusions and Critiques #
I would guess from this analysis that the sex differences have an effect on the methylation status of IDH, and possibly also TET2, but as there are so few probes it is difficult to say this with any certainty for TET2. This low number also suggests to me that we do not see GW hypermethylation of TET2? 
  It is also indicated that the %purity of obtained samples has a larger effect on methlyation status of both genes than either sex differences or mutation. Hence, it would be useful to do the mutational analysis with a threshold/ cut off for % purity of samples for both IDH and TET2 (looking at the design 3 table they don't always overlap very well). 


This analysis also suggests that the authors assumption of mutual exclusivity between IDH and TET2 mutations is not seen in their sample so further investigation into this group (4 from heatmap), would be beneficial.
  The authors are suggesting that LOF TET2 mutations are associated with the same epigenetic defects as IDH mutations, however using Bayesian analyses (fit 1) the mutant forms seemed to have wildly different methylation status (6513 for IDH and 6 for TET2), so I would like to look into this a bit further and see what they were getting at here! 


In sum, there are more probes for IDH than TET2, on the order of 1000x more probes no matter which design (controlling for sex, purity, etc.)
If this is the case, then the conclusions of this paper (that both mutations are associated with a hypermethylation pattern) should be called into question.
Instead, it appears that the methylation patterns are different between TET2 vs IDH mutated tumors.

