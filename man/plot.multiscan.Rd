%                              ~*~ Mode: Rd ~*~
%###############################################################################
%##                            plot.multiscan.Rd
%###############################################################################
%## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
%## Created On      : 2008-06-18 14:39
%## Last Modified By: Mizanur Khondoker
%## Last Modified On: 2008-06-19 22:55
%## Update Count    : 4
%###############################################################################
%##
%## Copyright (C) 2008 Mizanur Khondoker
%## Distributed under GPL 2 or later
%###############################################################################

\name{plot.multiscan}
\alias{plot.multiscan}

\title{Plot fitted model and standardised residuals for a multiscan object}
\description{
  This function provides  plots for fitted model and standardised
  residual for \code{\link{multiscan}} fitted objects.
}
\usage{
\method{plot}{multiscan}(x, residual=FALSE,\dots)
}

\arguments{
  \item{x}{a \code{multiscan} fitted object.}
  \item{residual}{Logical. Should residuals be plotted instead of the
    fitted model?. Defaults to \code{FALSE}.}
   \item{\dots}{Further arguments passed to the \code{plot} function.}
}
\details{
  Fitted model on the input data, after rescaling by the corresponding
  scanning effects, is plotted against the estimated gene expressions.
  Standardised residuals for each scan of data are plotted against the
  rank of estimated gene expressions.} 

\value{
   Returns either one plot of the fitted model (\code{residual=FALSE}) on the input data or
   \code{m} plots of the residuals (\code{residual=TRUE}) corresponding to each scan of data.}

\references{ 

 Khondoker, M. R., Glasbey, C. A. and Worton, B. J. (2006).
  Statistical estimation of gene expression using multiple laser scans of microarrays.
  \emph{Bioinformatics} \bold{22}, 215--219.}


\seealso{ \code{\link{multiscan}} }

\examples{

data(murine)

fit<-multiscan(murine)

## plot the fitted model

plot(fit)

## plot the residuals

op<-par(mfrow=c(2,2))
plot(fit, residual=TRUE)
par(op)
}
\keyword{regression} 
