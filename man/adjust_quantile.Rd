\name{adjust_quantile}
\alias{adjust_quantile}

\title{
A function to compute the quantiles of pmax using either approximation
or exact method
%%  ~~function to do ... ~~
}
\description{
A function to generate the quantiles of pmax 
(the maximum of the two p-values using the mixture null distribution)
}
\usage{
adjust_quantile(alpha00, alpha01, alpha10, alpha1, alpha2, 
input_pvalues, method = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha00}{
A numeric number represents the proportion of null \eqn{H_{00}}
}
  \item{alpha01}{
A numeric number represents the proportion of null \eqn{H_{01}}
}
  \item{alpha10}{
A numeric number represents the proportion of null \eqn{H_{10}}
}
  \item{alpha1}{
A numeric number represents the proportion of null alpha==0 
(association between exposure and mediator)}
  \item{alpha2}{
A numeric number represents the proportion of null beta==0 
(association between mediator and outcome adjusted for exposure)}
  \item{input_pvalues}{
A matrix contains two columns of p-values for candidate mediators.
Column 1 is the p-value of testing if a exposure is associated with 
the mediator (alpha!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome adjusted for exposure(beta!=0)  
}
  \item{method}{
Use the option to choose from two methods. 
method=0: approximation method without estimating the CDFs;
method=1:estimate the CDFs nonparametrically}
}
\details{
A function used to generate the quantiles for the observed 
p-max values using the mixture null distribution.
}
\value{
A vector contains the p-values based on generated null distribution

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
data(snp_input)
input_pvalues=snp_input
#To save time for illustration, we use 10 percent of rows
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),
size=ceiling(nrow(input_pvalues)/10)),]

nullprop <- nullestimation(input_pvalues,lambda=0.5)
pnull <- adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,method=0)
}

\keyword{composite null}
\keyword{joint significance}
\keyword{mediation}