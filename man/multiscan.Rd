%                              ~*~ Mode: Rd ~*~
%###############################################################################
%##                              multiscan.Rd
%###############################################################################
%## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
%## Created On      : 2008-04-22 14:39
%## Last Modified By: Mizanur Khondoker
%## Last Modified On: 2008-06-07 16:52
%## Update Count    : 8
%###############################################################################
%##
%## Copyright (C) 2008 Mizanur Khondoker
%## Distributed under GPL 2 or later
%###############################################################################

\name{multiscan}
\alias{multiscan}
\title{Combining multiple laser scans of microarray data}
\concept{minimization}
\concept{maximization}
\description{
  Estimates gene expressions from multiple laser scans of microarrays using non-linear 
  functional regression model with additive plus multiplicative errors.
}
\usage{
multiscan(data, initial = NULL, na.rm = TRUE, verbose = FALSE, control = list())
}

\arguments{
 \item{data}{A numeric matrix or data frame containing the intensity data 
 of a single microarray scanned at multiple (two or more) scanner settings. For
dual channel arrays, \code{multiscan} should be run on each channel of data separately.
 Number of rows (\code{n}) is equal to the number of spots/probes on the array, 
 and the number of columns (\code{m}) equals the number of scans. Columns will be
 arranged in order of scanner's sensitivity before fitting the model. Replicated 
 probes on the array  are treated as individual spots.}
 
 \item{initial}{A vector of length \code{m+2} to be used as initial values for 
 the scanning effects $(\beta_2,\cdots, \beta_m)$ and scale 
 $(\sigma_1, \sigma_2, \nu)$  parameters. If it is \code{NULL} (default), the 
 initial values are determined from the \code{data}.}
 
 \item{na.rm}{Logical. Should missing values be removed? Defaults to \code{TRUE}.}
 
 \item{verbose}{Logical. If \code{TRUE}, some intermediate results are printed 
 as the iteration proceeds.}
 
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 
 
}
\details{
  The function implements the method of Khondoker \emph{et. al.} (2006) 
  for combining multiple laser scans of microarrays. This function is 
  computationally slow and memory-intensive. That is due to the nested
  iteration loops of the numerical optimisation of the likelihood
  function involving a large number ($n+m+2$) of parameters. The 
  optimisation uses an alternating algorithm with Nelder-Mead simplex
  (Nelder and Mead, 1965) in the inner loops. The function \code{multiscan}
  directly uses the C function \code{nmmin}, the internal code used in the 
  general-purpose optimisation tool \code{\link{optim}}, for implementing Nelder-Mead
  simplex. For large arrays with many tens of thousands of probes, it is recommended 
  to consider a random subset (10-15,000 rows) of the data matrix to fit an initial model, 
  and then use the estimated $\beta$ and scale parameters as initial values for fitting
  the model on the full data set.  
  
  The \code{control} is a list of arguments. The users can change/supply any of the
  following components:
  \describe{
    \item{\code{trace}}{Indicator (\code{0} or \code{1}) of tracing information of 
    Nelder-Mead algorithm. If \code{1}, tracing information on the progress of the 
    optimisation is produced. Because Nelder-Mead may be callled thousands of times 
    during the estimation process, setting \code{trace = 1} will print too much 
    information very rapidly, which may not be useful. Defaults to \code{0}.}
    
     \item{\code{gmaxit}}{The maximum number of global iterations. Defaults to \code{150}.}  
    
     \item{\code{maxit}}{The maximum number of Nelder-Mead iterations. 
     Defaults to \code{5000}.}  
    
    \item{\code{reltol}}{Relative convergence tolerance of Nelder-Mead.  
    The algorithm stops if it is unable to reduce the value by a factor of
      \code{reltol * (abs(val) + reltol)} at a step.  Defaults to \code{1e-5}.}
      
    \item{\code{globaltol}}{Convergence tolerance of the outer (alternating) iteration.
    The estimation process converges if the gain in loglikelihood from one complete cycle
    of the outer iteration is less than \code{globaltol}. Defaults to \code{1e-10}.} 
          
    
    \item{\code{alpha}, \code{beta}, \code{gamma}}{Scaling parameters
      for the Nelder-Mead method. \code{alpha} is the reflection
      factor (default 1.0), \code{beta} the contraction factor (0.5) and
      \code{gamma} the expansion factor (2.0).}
    }

 
}
\value{
   Returns an object of class \code{multiscan} with components
  \item{call}{The call of the \code{multiscan} function.}

  \item{beta}{A vector of length \code{m} containing the maximum likelihood estimates of 
  the scanning effects, the first component fixed at \code{1}.}
  
  \item{scale}{A vector of length \code{3} containing the maximum likelihood estimates of 
  the scale parameters $\sigma_1, \sigma_2, \mbox{and} \nu$.}
  
  \item{mu}{A vector of length \code{n} containing the estimated gene expressions.}
  \item{data}{A matrix of the input data with columns rearranged in order of scanner's sensitivity.}
 \item{fitted}{A matrix of the fitted model on the \code{data}.}
 \item{sdres}{A matrix of the standardised residuals.}
  \item{outerit}{Number of  global iterations completed.}
 
 \item{gconv, conv, convmu}{Integer convergence codes. 
 
     \describe{
      \item{\code{gconv}}{Indicator of global convergence. \code{0} indicates successful convergence,
      \code{1} indicates premature termination.}
      \item{\code{conv}}{Convergence codes for Nelder--Mead simplex in the last global iteration while 
            updating scanning effects and scale parameters. \code{0} for successful convergence,
	    \code{1} indicates that the iteration limit \code{maxit} had been reached, \code{10} 
	     indicates degeneracy of the Nelder--Mead simplex.}
	\item{\code{convmu}}{Convergence codes for Nelder--Mead simplex in the last global iteration while 
            updating the gene expression parameters. This is an integer vector of length \code{n} where
	    each component takes the value \code{0}, \code{1}, or \code{10} depending on whether the 
	    Nelder--Mead simplex successfully converged, reached iteration limit \code{maxit} or produced
	    degeneracy respectively while updating the corresponding gene expression parameter.}
	      
	}
  }
 
  \item{\code{outerit}}{Number of  global iterations completed.}
 
 \item{\code{loglf}}{Value of the loglikelihood function at convergence (\code{gconv=0}). \code{NA}
 if not converged (\code{gconv=1}).}	
 
}

\references{
  Khondoker, M. R., Glasbey, C. A. and Worton, B. J. (2006).
  Statistical estimation of gene expression using multiple laser scans of microarrays.
  \emph{Bioinformatics}, \bold{22}, 215--219.

  Nelder, J. A. and Mead, R. (1965). 
  A simplex algorithm for function minimization. 
  \emph{Computer Journal} \bold{7}, 308--313.
 }
 
 \seealso{
 A web interface, created by David Nutter of Biomathematics & Statistics Scotland (BioSS), based on the original 
 Fortran code is available at http://www.bioss.ac.uk/ktshowcase/create.cgi. Results from the web interface may 
 not be exactly identical to that of \code{multiscan} as it uses a different (non-free IMSL routine) implementation 
 of Nelder--Mead simplex. 
 }
 

\examples{

## load the multiscan library 
library(multiscan)

## load the murine data set included in multiscan package
data(murine)
murine[1:10,] ## see first few rows of data

## fit the model on murine data with default options
fit<-multiscan(murine)
fit

## plot the fitted model
plot(fit)

## get the estimated gene expressions
gene.exprs<-fit$mu

## see more details as iteration progresses

fit1<-multiscan(murine, verbose = TRUE)
fit1

}
\keyword{nonlinear}
\keyword{optimize}
