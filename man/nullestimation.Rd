\name{nullestimation}
\alias{nullestimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the proportion of three component nulls
}
\description{
A function to estimate the proportion of three component
nulls. We developed a JS (joint significance)-mixture method based the 
estimated nulls which can provide much more accurate control of the 
family-wise error rate (FWER) and the false discovery rate (FDR) compared to
the JS-uniform method.
}
\usage{
nullestimation(input_pvalues, lambda = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_pvalues}{
A matrix contains two columns of p-values for candidate mediators.
Column 1 is the p-value of testing if an exposure is associated with 
the mediator (alpha!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome adjusted for the exposure (beta!=0) 

}
  \item{lambda}{
A tuning parameter between 0 and 1, the default value is 0.5.
}
}
\details{
A function used to estimate  the proportions of the three types 
of component null hypotheses:

\eqn{H_{00}}: alpha==0 and beta==0

\eqn{H_{01}}: alpha==0 and beta!=0

\eqn{H_{10}}: alpha!=0 and beta==0
}
\value{
A list contains five elements.
\item{alpha00}{A numeric number represents the proportion of null \eqn{H_{00}}}
\item{alpha01}{A numeric number represents the proportion of null \eqn{H_{01}}}
\item{alpha10}{A numeric number represents the proportion of null \eqn{H_{10}}}
\item{alpha1}{A numeric number represents the proportion of null beta==0}
\item{alpha2}{A numeric number represents the proportion of null alpha==0}
}

\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2019, submitted.  
}
\author{
James Y. Dai and X. Wang
}

\examples{
data(exampleinputdata)
input_pvalues=snp_input
#To save time for illustration, we use 10 percent of rows
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),
size=ceiling(nrow(input_pvalues)/10)),]

nullprop <- nullestimation(input_pvalues,lambda=0.5)
}

\keyword{composite null, joint significance, mediation}
