## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(HDMT)

## ----data---------------------------------------------------------------------
data(snp_input)
dim(snp_input)
data(exercise_input)
dim(exercise_input)

## ----input pvalue-------------------------------------------------------------
input_pvalues <- snp_input
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),size=ceiling(nrow(input_pvalues)/10)),]

## -----------------------------------------------------------------------------
nullprop <- null_estimation(input_pvalues)
nullprop

## -----------------------------------------------------------------------------
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)

## ----fig.height=5,fig.width=5-------------------------------------------------
pmax <- apply(input_pvalues,1,max)
correct_qqplot(pmax, pnull1)

## ----fig.height=5,fig.width=5-------------------------------------------------
efdr <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)

plot(pmax[order(pmax)],efdr[order(pmax)],type="l",ylim=c(0,1),xlab="p-max",ylab="Estimated FDR")

## -----------------------------------------------------------------------------
input_pvalues <- exercise_input
#To save time, we use 10% of rows
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),size=ceiling(nrow(input_pvalues)/10)),]

## -----------------------------------------------------------------------------
nullprop <- null_estimation(input_pvalues)

## -----------------------------------------------------------------------------
pnull<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)

## -----------------------------------------------------------------------------
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)

## ----fig.height=5,fig.width=5-------------------------------------------------
pmax <- apply(input_pvalues,1,max)
correct_qqplot(pmax, pnull)

## ----fig.height=5,fig.width=5-------------------------------------------------
correct_qqplot(pmax, pnull1)

## ----sessionInfo, include=TRUE, echo=TRUE, results='markup'-------------------
sessionInfo()

