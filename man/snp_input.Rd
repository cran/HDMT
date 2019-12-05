\name{snp_input}
\alias{snp_input}
\docType{data}
\title{
An example dataset to demonstrate the usage of 'HDMT'
}
\description{
The dataset was included to show the mediation role of DNA methylation 
in the study of genetic regulation of gene expression in primary prostate 
cancer (PCa) samples from The Cancer Genome Atlas (TCGA) with risk SNPs 
as the exposure.

The dataset is a matrix containing two columns of p-values for 
candidate mediators. Column 1 is the p-value of testing if a mediator 
is associated with the outcome conditional on the exposure (beta!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome (alpha!=0). 
}
\usage{data("snp_input")}
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
data(snp_input)
dim(snp_input)

}
\keyword{dataset}
