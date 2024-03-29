% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{computeUnivariateTernaryMatrix}
\alias{computeUnivariateTernaryMatrix}
\title{Compute the ternary matrix with digitized divergence coding}
\usage{
computeUnivariateTernaryMatrix(Mat, Ranges)
}
\arguments{
\item{Mat}{Matrix of data to be digitized (usually in quantiles), with each column corresponding to a sample and each 
row corresponding to a feature}

\item{Ranges}{A data frame containing the baseline range of each features; this corresponds to the 'Ranges' element of 
the output of findUnivariateGammaWithSupport() and computeUnivariateSupport()}
}
\value{
A matrix with the same dimensions as Mat containing the ternary form data.
}
\description{
Function for obtaining the ternary form for a matrix of data given a baseline range
}
\examples{
data(divergence)
S = computeUnivariateSupport(computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group=="NORMAL"]))
T = computeUnivariateTernaryMatrix(Mat=computeQuantileMatrix(breastTCGA_Mat[, breastTCGA_Group!="NORMAL", ]), Ranges=S$Ranges)
}
\keyword{digitization}
\keyword{ternary}
