##                             ~*~ Mode: R ~*~
################################################################################
##                                  zzz.R
################################################################################
## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
## Created On      : 2008-05-30 14:39
## Last Modified By: Mizanur Khondoker
## Last Modified On: 2008-05-31 12:31
## Update Count    : 9
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

#.onLoad<-function(lib,pkg){
#    require("Biobase",  quietly=TRUE) || stop("'Biobase' package not found")
#}


.onAttach<-function(lib, pkg)
{
    cat(paste("*** Loaded multiscan Version",
              packageDescription("multiscan")$Version), "*** \n")
    cat("Type 'vignette(\"multiscan\")' to view the package vignette.\n")

    ## to put vignette to the R menu bar (for Windows/GUI users)
    
    if(interactive() && .Platform$OS.type=="windows" && .Platform$GUI=="Rgui"){
        addVigs2WinMenu("multiscan")
    }
}
