#' Fit linear models and/or contrasts to data from Figueroe et al. 2010. 
#' 
#' @param   eset    an ExpressionSet of HELP data (e.g. data(DNAme))
#' @param   design  a formula for a linear model to fit to the data with limma
#' @param   ...     any other parameters to pass to lmFit (e.g. method="robust")
#'
#' @return          results from eBayes(lmFit(eset, design, ...))
#' 
#' @examples
#' 
#'  data(DNAme, package="WorldsSimplestCodeReview")
#'  covariates <- pData(DNAme) # pData == "pheno-data", a very old function
#'  library(limma) 
#'
#'  # let's fit a model with the specified bone marrow blast % (but beware...)
#'  design <- with(covariates, model.matrix( ~ IDH:purity + TET2:purity + male))
#'  fitWithPurityAndSex <- TET2(DNAme, design, method="robust")
#' 
#'  # if you look at DNAme$purity, it doesn't make a lot of sense
#'  # so let's compare against normal bone marrows (NBM) instead:
#'  DNAme$NBM <- grepl("NBM", DNAme$title) 
#'  DNAme$AML <- !DNAme$NBM # the non-NBMs are all AMLs 
#'  design2 <- with(pData(DNAme), model.matrix( ~ IDH:AML + TET2:AML + AML))
#'  fitAgainstNBMs <- TET2(DNAme, design2) # AML-specific hyper now in coef 4!
#'  colnames(design2) # for reference below
#' 
#'  message("IDH vs. non-mutant AML @ p_adj < 0.05:")
#'  nrow(topTable(fitAgainstNBMs, coef=2, p.val=0.05, n=Inf)))
#' 
#'  message("TET2 vs. non-mutant AML @ p_adj < 0.05:")
#'  nrow(topTable(fitAgainstNBMs, coef=3, p.val=0.05, n=Inf))
#' 
#'  message("AML vs. normal bone marrow @ p_adj < 0.05:")
#'  nrow(topTable(fitAgainstNBMs, coef=4, p.val=0.05, n=Inf))
#' 
#' @import          limma
#' 
#' @export
TET2 <- function(eset, design, ...) {
 
  message("Fitting...")
  fit0 <- lmFit(eset, design, ...) 
  
  message("Applying empirical Bayes shrinkage...")
  fit1 <- try(eBayes(fit0))

  # if shrinkage worked... 
  if (!inherits(fit1, "try-error")) {
    
    # and if there's a TET2 coefficient... 
    if (grepl("TET2", colnames(design))) {

      # enumerate the hits for it from the fit
      TET2coef <- grep("TET2", colnames(design))[1]
      message("Enumerating hits for ", 
              colnames(design)[TET2coef], 
              " at p_adj < 0.05 ...")

      # and stuff them into an attribute of the result for later
      tt <- topTable(fit1, coef=TET2coef, p.value=0.05, number=Inf)
      attr(fit1, "TET2toptable") <- tt

      # then count them off 
      message(nrow(tt), " probes significant for ", 
              colnames(design)[TET2coef], " at p_adj < 0.05")

    }

    # regardless, return the shrunken fit:
    message("Returning empirical Bayes shrunken model fit results.") 
    return(fit1) 

  } else { 

    warning("eBayes shrinkage failed, returning raw lmFit results.")
    return(fit0) 

  }

}
