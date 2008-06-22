##                             ~*~ Mode: R ~*~
################################################################################
##                             plot.multiscan.R
################################################################################
## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
## Created On      : 2008-03-22 14:30
## Last Modified By: Mizanur Khondoker
## Last Modified On: 2008-06-21 20:10
## Update Count    : 8
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

plot.multiscan<-function(x, residual=FALSE,...)
{
    if(!inherits(x,"multiscan"))
        stop("This is not a \'multiscan\' object")
    m<-dim(x$data)[2]
   
    ## Plot the fitted model
    ## -----------------------------------------------------------------------

    if (!residual){
        ordmu<-order(x$mu)
        colmax<-apply(x$data,2,max)
        ymax<-max(colmax/x$beta)
        xlab<-expression(paste("Estimated gene expression ",(hat(mu)[i]),sep=""))
        ylab<-expression(paste("Rescaled intensities ",(y[ij]/hat(beta)[j]),sep=""))
        par(mar=c(5, 5, 4, 2)+ 0.1)
        plot(x$mu,x$data[,1]/x$beta[1],ylim=c(0,ymax),type="n", xlab=xlab,ylab=ylab,...)
    
        for (j in 1:m){
            points(x$mu,x$data[,j]/x$beta[j], pch=j,col=j)
            lines(x$mu[ordmu],x$fitted[,j][ordmu]/x$beta[j],lty=2,col=j)
        }
        legend(0,0.9*ymax, colnames(x$data),pch=c(1:m),col=c(1:m),lty=rep(2,m))
    }

    ## Residual plot
    ## -----------------------------------------------------------------------
    if (residual)
    {
        sdres<-x$sdres
        mu<-x$mu
        for (i in 1:m){
            if(i==2)op<-par(ask = prod(par("mfcol"))<m && dev.interactive())
            plot(rank(mu), sdres[,i], xlab="Rank of estimated gene expression", 
                 ylab="Standardised residuals",main=colnames(x$data)[i], cex=0.8,...)
            abline(h=0, col="green")
            abline(h=-12.71, col="blue", lty=2)
            abline(h=12.71, col="blue", lty=2)
        }
        par(op)
    }
}

## ---------------------------------------------------------------------------
## End:     plot.multiscan.R
##----------------------------------------------------------------------------  
  
