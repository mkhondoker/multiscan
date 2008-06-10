/*                             ~*~ Mode:C ~*~
################################################################################
##                              multiscan.c
################################################################################
## Author          : Mizanur Khondoker, Chris Glasbey, Bruce Worton
## Created On      : 2008-03-21 21:18
## Last Modified By: Mizanur Khondoker
## Last Modified On: 2008-05-06 08:53
## Update Count    : 14
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
*/

#include <R.h>
#include <Rmath.h>        
#include <R_ext/Applic.h>  
                            
double 	*mem_vec(int nc);
double 	**mem_mat(int nr, int nc);

void    vec2vec(double *invec, double *outvec, int nc);
void 	vec2mat(double *vec, double **mat, int nr, int nc);
void  	mat2vec(double **mat, double *vec, int nr, int nc);

#define T 65535.0  /* the censoring threshold */ 


typedef struct{

    double *yi;
    double *b;
    double *s1;
    double *s2;
    double nu;
    double t;
    int m;
}data1;


typedef struct{
    
    double **y;
    double *mu;
    double *scale;
    double t;
    int n;
    int m;
}data;

double loglik1(int np1, double *par1, void *ex1)
{
    int j,lowerTail=1,giveLog=0;
    double fnvalue1=0.0;
    double mean=0.0,sd=1.0;
    double mui;
    double cm1, neglf1;
 
    data1 *fixed1=ex1;
    double *yi=fixed1->yi;
    double *b=fixed1->b;
    double *s1=fixed1->s1;
    double *s2=fixed1->s2;
    double nu=fixed1 ->nu;
    double t=fixed1->t;
    int m=fixed1->m;
    
    mui=par1[0];
    for (j=0;j<m;j++){
	cm1=t+(mui*b[j]-t)*pnorm((t-mui*b[j])/(fabs(mui*b[j]*nu)),mean,sd,lowerTail,giveLog)-
	    (fabs(mui*b[j]*nu))*dnorm((t-mui*b[j])/(fabs(mui*b[j]*nu)),mean,sd,giveLog);
        neglf1=0.5*log(pow(s1[j],2)+ pow(s2[j]*mui,2))+
	    log (1.0+ (pow(yi[j]-cm1,2))/(pow(s1[j],2)+ pow(s2[j]*mui,2)));
        fnvalue1+=neglf1;         
    }
    return fnvalue1;
}

double loglik(int np, double *par, void *ex)
{
    int i,j,k,lowerTail=1,giveLog=0;
    double fnvalue=0.0, sigma1, sigma2, nu;
    double mean=0.0,sd=1.0;
    double *para, *b, *s1, *s2;
    double cm, neglf;
 
 

    data *fixed=ex;
    double **y=fixed->y;
    double *mu=fixed->mu;
    double *scale=fixed ->scale;
    double t=fixed->t;
    int n=fixed->n;
    int m=fixed->m;

    para=mem_vec(np);
    b=mem_vec(m);
    s1=mem_vec(m);
    s2=mem_vec(m);
 
    for (k=0; k<np; k++){
	para[k]=par[k]*scale[k];
    }

    sigma1=para[m-1];
    sigma2=para[m];
    nu=para[m+1]; 

    for (j=0;j<m;j++){
	if (j==0){
	    b[j]=1.0;
	}
	else{
	    b[j]=para[j-1];
	}
	s1[j]=sigma1*b[j];
	s2[j]=sigma2*b[j];
    }
    
    for (i=0;i<n;i++){
	for (j=0;j<m;j++){
	    cm=t+(mu[i]*b[j]-t)*pnorm((t-mu[i]*b[j])/(fabs(mu[i]*b[j]*nu)),mean,sd,lowerTail,giveLog)-
		(fabs(mu[i]*b[j]*nu))*dnorm((t-mu[i]*b[j])/(fabs(mu[i]*b[j]*nu)),mean,sd,giveLog);
	    neglf=0.5*log(pow(s1[j],2)+ pow(s2[j]*mu[i],2))+
		log (1.0+ (pow(y[i][j]-cm,2))/(pow(s1[j],2)+ pow(s2[j]*mu[i],2)));
	    fnvalue+=neglf;         
	}
    }
 
    return fnvalue;
}


void multiscan (int *n, int *m, int *np, double *y, double *mu, double *initial, double *scale, 
		double *est,double *fitted,double *sdres, int *fail, double *reltol,double *globaltol, 
		double *alpha,	double *beta, double *gamma, int *trace,int *verbose,int *gmaxit,
		int *maxit, int *outerit,int *gconv,int *convmu, double *loglf)

{ 
    int i,j,k,np1=1,bigint=300000,jmin,failmu,fncount;
    double dk,fmin,sdij;
 
    double initlf,oldlf,newlf,gainlf, sigma1, sigma2,locnu,minlfm,abstol=-1.0e35;
    double *par,*par1,*estimate,*muest,*loc,*mu1,*mfailmu;
    double **mmu,**mmu1,**mlfm,**mfitted,**msdres;
  

    par=mem_vec(*np);
    par1=mem_vec(np1);
    estimate=mem_vec(*np);
    muest=mem_vec(np1);
    mu1=mem_vec(*n);
    loc=mem_vec(bigint);
    mmu=mem_mat(*n,*m);
    mmu1=mem_mat(*n,*m);
    mlfm=mem_mat(*n,*m);
    mfailmu=mem_vec(*m);
    mfitted=mem_mat(*n,*m);
    msdres=mem_mat(*n,*m);
 
    data1 ex1;
    ex1.yi=mem_vec(*m);
    ex1.b=mem_vec(*m);
    ex1.s1=mem_vec(*m);
    ex1.s2=mem_vec(*m);
    ex1.t=T;
    ex1.m=*m;
 
    data ex;
    ex.y=mem_mat(*n,*m);
    ex.mu=mem_vec(*n);
    ex.scale=mem_vec(*np);
    ex.t=T;
    ex.n=*n;
    ex.m=*m;

    vec2mat(y,ex.y,*n,*m);
    for (i=0;i<*n;i++){
	ex.mu[i]=mu[i];
	if(ex.mu[i]<=0.0){
	    ex.mu[i]=1.0e-2;
	}
    }
    for (k=0;k<*np;k++){
	ex.scale[k]=scale[k];
	par[k]=initial[k]/scale[k];
    }
  
  /* printmat(ex.y, *n,*m);*/
    Rprintf("\n");
    Rprintf("Estimating %d gene expressions using %d scans of data\n",*n,*m);
    if(*verbose){
	Rprintf("------------------------------------------------------------\n");
	Rprintf("Initial parameters:\n");
	Rprintf("------------------------------------------------------------\n");

	for (j=0;j<*m;j++){
	    if (j==0){
		ex1.b[j]=1.0;
	    }
	    else{
		ex1.b[j]=initial[j-1];
	    }
	}   
	
	sigma1=initial[*m-1];
	sigma2=initial[*m];
	ex1.nu=initial[*m+1];
	Rprintf("Scanning effects:\n");
	for (j=0;j<*m;j++){
	    Rprintf("beta[%d]=%14.8f\n",j+1,ex1.b[j]);
	}
	Rprintf("\n");
	Rprintf("Scale parameters:\n");
	Rprintf("sigma1=%14.8f\n",fabs(sigma1));
	Rprintf("sigma2=%14.8f\n",fabs(sigma2));
	Rprintf("    nu=%14.8f\n",fabs(ex1.nu)); 
    }
    Rprintf("------------------------------------------------------------\n");
    
    /*******************************************************************************
  Evaluate and print the negative log-likelihood at initial parameters
    ******************************************************************************/
    initlf=loglik(*np,par,&ex);
    if(!R_FINITE(initlf)){
	error("log-likelihood function can not be evaluated at initial parameters");
    }
    Rprintf("Log-likelihood at initial parameters: %12.4f\n",-initlf);
    Rprintf("------------------------------------------------------------\n");
    
    /*******************************************************************************
 Update the other parameters (except mu) with  mu initialised at the scan1 data.
    ********************************************************************************/
    

    gainlf=1.0e35; /* an arbitrary big number */
    *outerit=0;     /* count the number of global iteration */
    *gconv=1;
    *loglf=0.0;
    
    /* Start the loop for global iteration*/
    
    while (gainlf>*globaltol && *outerit<*gmaxit){ 
	if(*verbose || *trace!=0){
	    Rprintf("************************************************************\n");
	    Rprintf("Started iteration: %5d\n",*outerit+1);
	}
	oldlf=loglik(*np,par,&ex);
	if(*verbose || *trace!=0){
	    Rprintf("Updating the scanning effects and scale parameters...");
	    if(*trace!=0){
		Rprintf("\n");
	    }
	}
	nmmin(*np,par, est,&fmin, loglik,fail,abstol,*reltol, &ex,*alpha,*beta,
	      *gamma, *trace,&fncount, *maxit);
	for (j=0;j<*np;j++){
	    par[j]=est[j];
	    estimate[j]=par[j]*scale[j];

	}


	sigma1=estimate[*m-1];
	sigma2=estimate[*m];
	if (*verbose || *trace!=0){
	    if(*trace==0){
		Rprintf(" done.\n");
	    }
	    else{
		Rprintf("Finished updating scanning effects and scale parameters.\n");
	    }
	}
	
	/*******************************************************************************
        Update the gene expression parameters (mu) one at a time holding the other
        parameters (beta, sigma1, sigma2, and nu) fixed at the updated values. 
	********************************************************************************/ 
 
 
	ex1.nu=estimate[*m+1];
	for (j=0;j<*m;j++){
	    if (j==0){
		ex1.b[j]=1.0;
	    }
	    else{
		ex1.b[j]=estimate[j-1];
	    }
	    ex1.s1[j]=fabs(sigma1*ex1.b[j]);
	    ex1.s2[j]=fabs(sigma2*ex1.b[j]);
   
	}

	locnu=ex1.nu;
	if(locnu>0.1){
	    locnu=0.1;
	}
   
	for (k=0;k<bigint;k++){
	    dk=(double)k;
	    loc[k]=ex1.t+(dk-ex1.t)*pnorm((ex1.t-dk)/(fabs(dk*locnu)),0.0,1.0,1,0)-
		(fabs(dk*locnu))*dnorm((ex1.t-dk)/(fabs(dk*locnu)),0.0,1.0,0);
	    /* Rprintf("loc[%d]=%lf\n",k,loc[k]);*/
 }   

	for (i=0;i<*n;i++){
	    for (j=0;j<*m;j++){

		for (k=0;(k<=bigint && ex.y[i][j]-loc[k]>1.0e-8);k++){
		    ;
		}
		dk=(double)k;
		mmu[i][j]=dk/ex1.b[j]; 
		/* Rprintf("mmu[%d][%d]=%lf b[%d]=%lf\n",i,j,mmu[i][j],j,ex1.b[j]);*/  
	    }
	}

	if(*verbose || *trace!=0){
	    Rprintf("\nUpdating the gene expression parameters (mu)...");
	    if(*trace!=0){
		Rprintf("\n");
	    }
	} 
	for (i=0;i<*n;i++){
	    
	    
	    ex1.yi=ex.y[i];
	    for (j=0;j<*m;j++){
		par1[0]=mmu[i][j];
		
		if(*trace!=0){
		    Rprintf("Updating mu[%d] with starting value[%d] ...\n",i+1,j+1);
		}
   
		nmmin(np1,par1, muest,&fmin, loglik1, &failmu,abstol,*reltol, &ex1,*alpha,
		      *beta,*gamma, *trace,&fncount, *maxit);
		mmu1[i][j]=muest[0];
		mlfm[i][j]=fmin;
		mfailmu[j]=failmu;
		/* Rprintf("Initial mu:%lf Estimated mu:%lf fmin:%15.10lf\n",mmu[i][j],mmu1[i][j], mlfm[i][j]);*/
	    }
	    jmin=0;
	    minlfm=mlfm[i][jmin];
	    for (j=1;j<*m;j++){
		if (mlfm[i][j]<minlfm){
		    minlfm=mlfm[i][j];
		    jmin=j;
		}
	    }
	    mu1[i]=mmu1[i][jmin];
	    convmu[i]=mfailmu[jmin]; 
	    /* Rprintf("minmu:%lf minlfm:%15.10lf\n", mu1[i],minlfm);*/
	    if(*trace!=0){
		Rprintf("Finished updating mu[%d].\n",i+1);
	    }
    
	}
	if(*verbose || *trace!=0){
	    if(*trace==0){
		Rprintf(" done.\n");
	    }
	    else{
		Rprintf("Finished updating gene expression parameters.\n");
	    }
	}
	*outerit+=1; 

	/* Update mu in the structure ex */

	for (i=0;i<*n;i++){
	    ex.mu[i]=mu1[i]; 
	    if(ex.mu[i]<=0.0){
		ex.mu[i]=1.0e-2;
	    }
	}

	/* Evaluate the negative loglikelihood funnction with updated par and mu */

	newlf=loglik(*np,par,&ex);

  
	gainlf=fabs(oldlf-newlf);
	if (gainlf<=*globaltol){
	    *gconv=0;
	    *loglf=-newlf;
	}
	if(*verbose){
	    Rprintf("************************************************************\n");
	}
	Rprintf("End of global iteration: %d, Log-likelihood: %14.5f\n",*outerit,-newlf);
	if(*verbose){
	    
	    Rprintf("Estimates of the scanning effects:\n");
	    for (j=0;j<*m;j++){
		if(j==0){
		    Rprintf("beta[%d]=%14.8f (fixed)\n",j+1,ex1.b[j]);
		}
		else{
		    Rprintf("beta[%d]=%14.8f\n",j+1,ex1.b[j]);
		}
	    }
	    Rprintf("\n");
	    Rprintf("Estimates of the scale parameters:\n");
	    Rprintf("sigma1=%14.8f\n",fabs(sigma1));
	    Rprintf("sigma2=%14.8f\n",fabs(sigma2));
	    Rprintf("    nu=%14.8f\n",fabs(ex1.nu)); 
	    Rprintf("------------------------------------------------------------\n");
  }

  }  /* End of loop for global iteration */


    for (i=0;i<*n;i++){
	for (j=0;j<*m;j++){
	    mfitted[i][j]=ex1.t+(ex.mu[i]*ex1.b[j]-ex1.t)*pnorm((ex1.t-ex.mu[i]*ex1.b[j])/
		(fabs(ex.mu[i]*ex1.b[j]*ex1.nu)),0.0,1.0,1,0)- (fabs(ex.mu[i]*ex1.b[j]*ex1.nu))*
		dnorm((ex1.t-ex.mu[i]*ex1.b[j])/(fabs(ex.mu[i]*ex1.b[j]*ex1.nu)),0.0,1.0,0);
	    sdij=sqrt(pow(ex1.s1[j],2)+pow(ex.mu[i]*ex1.s2[j],2));
	    msdres[i][j]=(ex.y[i][j]-mfitted[i][j])/sdij;
     
	}

    }
    /* Copy 'mfitted' to 'fitted'; 'msdres' to 'sdres' */

    mat2vec(mfitted,fitted,*n,*m); 
    mat2vec(msdres,sdres,*n,*m); 

    vec2vec(ex.mu,mu,*n);
    vec2vec(estimate,est,*np);

} /* End of the main function "multiscan" */
 



double **mem_mat(int nr, int nc)

{
    int i;
    double **a;
    a=(double **)R_alloc( (unsigned) nr, sizeof(double *));
    if (a==NULL){
	error("Memory allocation failure 1 in mem_matrix()");
	return(NULL);
    }
    
    for (i=0; i<nr; i++){
	a[i]=(double *)R_alloc( (unsigned) nc, sizeof(double));
	if(a[i]==NULL){
	    error("Memory allocation failure 2 in mem_matrix()");
	    return(NULL);
	}
    }
    
    return a;
}


double *mem_vec(int nc)

{
    double *v;
    v=(double *)R_alloc( (unsigned) nc, sizeof(double));
    if (v==NULL){
	error("Memory allocation failure in mem_vec()");
	return(NULL);
    }
    return v;
}


void  vec2vec(double *invec, double *outvec, int nc)
    
{
    int i;
    for (i=0;i<nc;i++){
	outvec[i]=invec[i];
    }
    
}

void  vec2mat(double *vec, double **mat, int nr, int nc)

{
    int i,j;
    for (i=0;i<nr;i++){
	for (j=0;j<nc;j++){
	    mat[i][j]=vec[j+i*nc];
	}
    }

}

void  mat2vec(double **mat, double *vec, int nr, int nc)

{
    int i,j;
    for (i=0;i<nr;i++){
	for (j=0;j<nc;j++){
	    vec[j+i*nc]= mat[i][j];
	}
    }

}


/*
void  setmat0(double **mat, int nr, int nc)

{
  int i,j;
  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
     mat[i][j]=0.0;
    }
  }

}
*/

/*
void  setvec0(double *vec, int nc)

{
  int j;
     for (j=0;j<nc;j++){
     vec[j]=0.0;
     }

}
*/

/*
void printmat (double **mat, int nr, int nc)

{
  int i,j;
  for (i=0; i<nr;i++){
    Rprintf("[%d,] ",i);
    for (j=0;j<nc;j++){
      Rprintf("%14.4f  ", mat[i][j]);
    }
    Rprintf("\n");
}
}
*/
