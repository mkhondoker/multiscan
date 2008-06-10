##                             ~*~ Mode: R ~*~
################################################################################
##                              multiscan.R
################################################################################
## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
## Created On      : 2008-03-22 14:39
## Last Modified By: Mizanur Khondoker
## Last Modified On: 2008-05-30 22:08
## Update Count    : 15
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


multiscan<-function(data,initial=NULL,na.rm=TRUE,verbose=FALSE, control=list())
{
    if(!is.matrix(data) && !is.data.frame(data)){
        stop("data is not a matrix or data frame")
    }
    if(any(data<0)){
        stop("negative value(s) found in data : intensity can not be negative") 
    }
    if(na.rm){
        data <- na.omit(data)
    }
    else{
        na.fail(data)
    }
    data<-as.matrix(data)
    n<-dim(data)[1]
    m<-dim(data)[2]
    np<-m+2
    if(is.null(colnames(data))){
        colnames(data)<-paste("Column",1:m,sep="")
    }
    med<-apply(data,2,median)
    idx.scan<-order(med)
    if(!identical(idx.scan,c(1:m))){
        data<-data[,idx.scan]
    }
    mu<-data[,1]  
    if (m<2){
        stop("Number of scans must be at least 2")
    }
    if (is.null(initial)){
        beta<-(med[idx.scan])[2:m]/(med[idx.scan])[1]
        initial<-scale<-c(beta,10.0,0.005,0.1)
    }
    else{
        na.fail(initial)
        if(length(initial)!=np){
            stop("inappropriate length of initial: should be of length ",np)
        }
        if(any(initial<=0)){
            stop("all parameters should have positive values")
        }
        scale<-initial
    }
  
    est<-vector("numeric",np)
    fail<-99
    gconv<-99
    outerit<-0 
    loglf<-0 
    convmu<-array(99,n)
    fitted<-array(0,dim=c(n,m))
    sdres<-array(0,dim=c(n,m))    
    
    con <- list(trace =0,gmaxit=150,maxit=5000, reltol = 1.0e-5,
                globaltol=1e-10, alpha = 1,beta = 0.5, gamma = 2)
    con[(namc <- names(control))] <- control

    if(con$trace<0)
        stop("inappropriate value of 'trace'")
    if(con$maxit<=0)
        stop("unrealistic value of 'maxit'")
    
    if(con$gmaxit<=0)
        stop("unrealistic value of 'gmaxit'")
    
    cl<-match.call()
    
    fit<-.C("multiscan", 
            as.integer(n),
            as.integer(m),
            as.integer(np),
            as.double(t(data)),
            mu=as.double(mu),
            as.double(initial),
            as.double(scale),
            est=as.double(est),
            fitted=as.double(fitted),
            sdres=as.double(sdres),
            conv=as.integer(fail),
            as.double(con$reltol),
            as.double(con$globaltol),
            as.double(con$alpha),
            as.double(con$beta),
            as.double(con$gamma),
            as.integer(con$trace),
            as.integer(verbose),
            as.integer(con$gmaxit),	
            as.integer(con$maxit),
            outerit=as.integer(outerit),
            gconv=as.integer(gconv),
            convmu=as.integer(convmu),
            loglf=as.double(loglf),
            PACKAGE="multiscan"
            )
    out<-list()
    out$call<-cl
    estpars<-fit[["est"]]
    b<-c(1,estpars[1:(m-1)])
    out$beta<-b
    sigma1<-estpars[m]
    sigma2<-estpars[m+1]
    nu<-estpars[m+2]
    out$scale<-c(sigma1,sigma2,nu)
    out$mu<-fit[["mu"]]
    out$data<-data
    out$fitted<-matrix(fit[["fitted"]],ncol=m,byrow=T)
    out$sdres<-matrix(fit[["sdres"]],ncol=m,byrow=T)
    out$outerit<-fit[["outerit"]]
    out$gconv<-fit[["gconv"]]  
    out$convmu<-fit[["convmu"]]
    out$conv<-fit[["conv"]]
    out$loglf<-fit[["loglf"]]
    
    if(out$gconv==1){
        out$loglf<-NA
        warning("global iteration has non-zero exit status")
        if(out$outerit>=con$gmaxit){
            warning("global iteration reached 'gmaxit' without convergence")
        }
    }
 
    if(out$conv!=0){
        if(out$conv==1){
            warning("Nelder-Mead reached 'maxit' while optimising 'loglik'")
        }
        else
        {
            warning("degeneracy in Nelder-Mead simplex while optimizing 'loglik'")
        }
    }
    
    if(sum(out$convmu)!=0){
        if(any(out$convmu==1)){
            warning("Nelder-Mead reached 'maxit' while optimising 'loglik1'")
        }
        if(any(convmu==10)) {
            warning("degeneracy in Nelder-Mead simplex while optimizing 'loglik1'")
        }
    }
    class(out)<-"multiscan"
    out
}
