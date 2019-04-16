\name{exampleinputdata}
\alias{exampleinputdata}
\docType{data}
\title{
Example datasets to demonstrate the usage of 'HDMT'
}
\description{
Two datasets were included to show the mediation role of DNA methylation 
in two studies:
snp_input: genetic regulation of gene expression in primary prostate cancer 
(PCa) samples from The Cancer Genome Atlas (TCGA) with risk SNPs 
as the exposure.
exercise_input: regulation of prostate cancer progression in a Seattle-based
cohort of patients diagnosed with clinically localized PCa with excercise 
as the exposure. Due to space limit, a subset (10 percent) of the 
full dataset used in the paper is included here.
Each dataset is a matrix containing two columns of p-values for 
candidate mediators. Column 1 is the p-value of testing if a mediator 
is associated with the outcome conditional on the exposure (beta!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome (alpha!=0). 
}
\usage{data("exampleinputdata")}
\format{
  The format of snp_input is:
 num [1:69602, 1:2] 0.106 0.999 0.101 0.173 0.89 ...
}

\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2019, submitted. 
}
\examples{
data(exampleinputdata)
dim(snp_input)
dim(exercise_input)
}
\keyword{dataset}
