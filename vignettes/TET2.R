## ----message=FALSE, loadpkgs, eval=FALSE-----------------------------------------------------------------
#> install.packages("remotes")
#> install.packages("BiocManager")
#> library(BiocManager)
#> if (!require("GEOquery")) {
#>   BiocManager::install("GEOquery")
#>   library(GEOquery)
#> }
#> if(!require("limma")) {
#>   BiocManager::install("limma")
#>   library(limma)
#> }
#> #Kate told me I needed this and then I started Rstudio again and now it
#> #seems to work.
#> BiocManager::install("VanAndelInstitute/WorldsSimplestCodeReview")
#> library(tidyverse)
#> library(knitr)
#> knitr::opts_chunk$set(echo = TRUE)
#> knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
#> if (!requireNamespace("BiocManager", quietly = TRUE))
#>     install.packages("BiocManager")
#> 
#> library(devtools)
#> load_all("./")


## ---- tangle, eval = FALSE, message = FALSE, echo = TRUE-------------------------------------------------
#> knitr::knit("TET2.Rmd", tangle = TRUE)
#> [1] "TET2.R"


## ---- fetchGEO-------------------------------------------------------------------------------------------

library(limma)
library(GEOquery)
if (!exists("DNAme")) data(DNAme)

if (FALSE) { # this takes about 5 minutes:

  # needed to fetch data
  library(GEOquery) 
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
# Features  Samples
#    25626      394



## ---- heatmap, eval=TRUE---------------------------------------------------------------------------------

# always plot your data
library(ComplexHeatmap)
mutations <- t(as.matrix(pData(DNAme)[, c("TET2", "IDH")]))
Heatmap(mutations, col=c("lightgray","darkred"), name="mutant", column_km=4,
        column_names_gp = gpar(fontsize = 7))



## ---- The OddBall----------------------------------------------------------------------------------------
library(tidyverse)
# one patient is the odd-ball here
as_tibble(DNAme$`idh1.idh2:ch1`) -> idh1_idh2
# since there is ch1 and 2, I compared both and they have the exact same information
# as_tibble(DNAme$`idh1.idh2:ch2`) -> idh1_idh2_next
# idh1_idh2 == idh1_idh2_next - returns TRUE
as_tibble(DNAme$`tet2:ch1`) -> tet
# as_tibble(DNAme$`tet2:ch2`) -> tet_2
# tet == tet_2 - returns TRUE
colnames(tet) <- c("TET")
colnames(idh1_idh2) <- c("IDH")
compiled <- cbind(tet, idh1_idh2)
View(compiled) # scrolled through and identified the patient/sample number that had the mutations in both TET2 and IDH


## ---- TET2_vs_IDH----------------------------------------------------------------------------------------

# model TET2 and IDH1/2 mutant related hypermethylation
# note: there are plenty of confounders (pb%, bm%, wbc) that could be included
library(limma) 

# simplest design
design1 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 ))
fit1 <- eBayes(lmFit(exprs(DNAme), design1))
(IDH_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("IDH", colnames(design1)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
# 6513 probes for IDH

(TET_diffmeth_probes_fit1 <- nrow(topTable(fit1, 
                                           coef=grep("TET2", colnames(design1)),
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
# 6 probes for TET2

# control for sex
design2 <- with(pData(DNAme), model.matrix( ~ IDH + TET2 + male ))
fit2 <- eBayes(lmFit(exprs(DNAme), design2))
(IDH_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                           coef=grep("IDH", colnames(design2)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
# 6651 probes for IDH 

(TET2_diffmeth_probes_fit2 <- nrow(topTable(fit2, 
                                            coef=grep("TET", colnames(design2)),
                                            p.value=0.05, # change if you like 
                                            number=Inf)))
# 7 probes for TET2

# control for blast count
design3 <- with(pData(DNAme), model.matrix( ~ IDH:purity + TET2:purity))
fit3 <- eBayes(lmFit(exprs(DNAme)[, as.integer(rownames(design3))], design3))

(IDH_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                           coef=grep("IDH", colnames(design3)), 
                                           p.value=0.05, # change if you like 
                                           number=Inf)))
# 7450 probes for IDH:purity

(TET2_diffmeth_probes_fit3 <- nrow(topTable(fit3, 
                                            coef=grep("TET", colnames(design3)),
                                            p.value=0.05, # change if you like 
                                            number=Inf)))
# 10 probes for TET2:purity


