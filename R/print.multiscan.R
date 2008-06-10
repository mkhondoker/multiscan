##                             ~*~ Mode: R ~*~
################################################################################
##                             print.multiscan.R
################################################################################
## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
## Created On      : 2008-03-22 14:30
## Last Modified By: Mizanur Khondoker
## Last Modified On: 2008-05-06 08:54
## Update Count    : 4
################################################################################
##
## Copyright (C) 2008 Mizanur Khondoker
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the “GNU General Public License” 
## as published by the Free Software Foundation; either version 2 
## of the License, or (at your option) any later version.
##  
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the “GNU General 
## Public License” for more details.
##   
## You may also obtain a copy of the “GNU General Public License” from 
## the Free Software Foundation by visiting their Web site or by writing to
##   
##
##  Free Software Foundation, Inc.
##  59 Temple Place - Suite 330
##  Boston, MA 02111-1307
##  USA
##
################################################################################

print.multiscan <- function (x,...) 
{  
    cat("\nCall: ", deparse(x$call),  "\n", fill=TRUE)
    cat("Scanning effects:\n")
    names(x$beta)<-colnames(x$data)
    print.default(format(x$beta),print.gap = 3,quote = FALSE)

    cat("\nScale parameters:\n")
    names(x$scale)<-c("Additive","Multiplicative","Nu")
    print.default(format(x$scale),print.gap = 3,quote = FALSE)

    cat("\nLog-likelihood at convergence:\n")
    names(x$loglf)<-"Loglf"
    print.default(format(x$loglf),quote = FALSE)
   
    invisible(x)
}

##---------------------------------------------------------------------------------
##    End: print.multiscan.R
##---------------------------------------------------------------------------------
