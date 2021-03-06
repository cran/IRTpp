% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploratory.R
\name{an.parallel}
\alias{an.parallel}
\title{Parallel Analysis}
\usage{
an.parallel(x = NA, iterations = 0, centile = 0, seed = 0, mat = NA,
  n = NA)
}
\arguments{
\item{x}{a matrix or a Dataframe that holds the test response data}

\item{iterations}{a number indicating the amount of iterations that 
representing the number of random data sets to be produced in the analysis.}

\item{centile}{a number between 1 and 99 indicating the centile used in estimating bias.}

\item{seed}{specifies that the random number is to be seeded with the supplied integer.}

\item{mat}{specifies that the procedure use the provided correlation matrix rather
than supplying a data matrix through x. The n argument must also be supplied when 
mat is used.}

\item{n}{the number of observations. Required when the correlation matrix is supplied 
with the mat option.}
}
\value{
Retained Components a scalar integer representing the number of components retained.

Adjusted eigenvalues a vector of the estimated eigenvalues adjusted.

Unadjusted eigenvalues a vector of the eigenvalues of the observed data from either
an unrotated principal component analysis.

Bias a vector of the estimated bias of the unadjusted eigenvalues
}
\description{
performs Horn's parallel analysis for a principal component.
}
\details{
Is a implementation of Horn's (1965) tecnique for evaluating the components retained
in a principle component analysis (PCA). This procedure is a adaptation of the
function paran of Package Paran.
}
\examples{
data <- simulateTest(model="2PL",items=10,individuals=1000) 
an.parallel(data$test, iterations = 100, centile = 99, seed = 12)
}
\references{
John L. Horn (1965). A rationale and test for the number of factors 
in factor analysis. Psychometrika, Volume 30, Number 2, Page 179.

Dinno A. 2009. Exploring the Sensitivity of Horn's Parallel Analysis to the 
Distributional Form of Simulated Data. Multivariate Behavioral Research. 44(3): 362-388
}

