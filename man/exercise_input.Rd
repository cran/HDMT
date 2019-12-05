\name{exercise_input}
\alias{exercise_input}
\docType{data}
\title{
An example dataset to demonstrate the usage of 'HDMT'
}
\description{
The dataset was included to show the mediation role of DNA methylation 
in the study of regulation of prostate cancer progression in a 
Seattle-based cohort of patients diagnosed with clinically localized 
PCa with excercise as the exposure. Due to space limit, 
a subset (10 percent) of the full dataset used in the paper is included here.

The dataset is a matrix containing two columns of p-values for 
candidate mediators. Column 1 is the p-value of testing if a mediator 
is associated with the outcome conditional on the exposure (beta!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome (alpha!=0). 
}
\usage{data("exercise_input")}
\format{
  The format of exercise_input is:
 num [1:47900, 1:2] 0.4966344 0.1048730 0.1005355 0.4946623 ...
}

\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2019, submitted. 
}
\examples{
data(exercise_input)
dim(exercise_input)

}
\keyword{dataset}
