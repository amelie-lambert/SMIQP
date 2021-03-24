/****************************************************************************
 //  This file is part of the src of "SMIQP",
 //  Copyright (C) 2019  Amélie Lambert
 //  All Rights Reserved.
 // 
 //     CEDRIC - CNAM 
 //     292 rue saint martin
 //     F-75141 Paris Cedex 03
 //     France
 //      
 //     amelie.lambert@cnam.fr    http://cedric.cnam.fr/~lamberta 
 // 
 //  This code is published under the Eclipse Public License.
 // Author: Amélie Lambert
 //****************************************************************************/

#include <stdio.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include<string.h>
#include <math.h>


#include"quad_prog.h"
#include "utilities.h"
#include <ilcplex/cplex.h>
#include"miqcr_cplex.h"

extern Q_MIQCP qqp;
extern long nodecount;
extern double EPS_LM;
extern double EPS_LOCAL_SOL;
extern double EPS_BETA;

/* parameters for CPLEX */
extern int MEMORY;
extern double ABS_GAP;
extern  double REL_GAP;
extern int VARSEL;
extern int TIME_LIMIT_CPLEX;
extern double * best_x;
extern double ub_bb;






int setqpproblemdata_miqcr_cplex (char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p, double **lb_p, double **ub_p, char ** ctype_p)
{
  char     *zprobname = NULL;     /* Problem name <= 16 characters */        
  double   *zobj = NULL;
  double   *zrhs = NULL;
  char     *zsense = NULL;
  int      *zmatbeg = NULL;
  int      *zmatcnt = NULL;
  int      *zmatind = NULL;
  double   *zmatval = NULL;
  double   *zlb = NULL;
  double   *zub = NULL;
  int      *zqmatbeg = NULL;
  int      *zqmatcnt = NULL;
  int      *zqmatind = NULL;
  double   *zqmatval = NULL;
  int      status = 0;
  char *zctype = NULL;
  

 
  int i,j,k,l,s,lim,lim_bis,lim_ter, s_bis;
   
  zprobname = alloc_string(16); 

  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
  zmatbeg   = alloc_vector(qqp->nb_col);   
  zmatcnt   = alloc_vector(qqp->nb_col); 
  
 
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i<= qqp->ind_max_x;i++)
    zmatbeg[i]= zmatbeg[i-1] + qqp->x_cont[i-1]+  qqp->m+qqp->p +qqp->mq+qqp->pq;
   
  /* y  */
  for(j=0;j<qqp->ind_max_int;j++)
    for(k=0;k<qqp->ind_max_x;k++)
      {
	if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)]))
	  if(j<=k)
	    {
	      zmatbeg[i]= zmatbeg[i-1]+ qqp->mq+qqp->pq+ 3 ;
	      i++;
	    }
	  else
	    {  
	      zmatbeg[i]= zmatbeg[i-1]+ qqp->mq+qqp->pq+ 2;
	      i++;
	    }
      }


  
  for(j=qqp->ind_max_int;j<qqp->ind_max_x;j++)
    for(k=0;k<qqp->ind_max_int;k++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)]))
	{
	  zmatbeg[i]= zmatbeg[i-1] + 1;
	  i++;
	}
  
  /*  z */
   for(l=0;l<qqp->nb_int;l++)
    {
      lim=log_base_2(qqp->u[l])+1;
      for(j=0;j<qqp->n;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(l,j,qqp->n)]))
	  for(k=0;k<lim;k++)
	    {
	      zmatbeg[i]= zmatbeg[i-1] + 5;
	      i++;
	    }
    }
  
  /*  t */
   for(l=0;l<qqp->nb_int;l++)
    {
      lim=log_base_2(qqp->u[l])+1;
      for(k=0;k<lim;k++)
	{
	  if (i < qqp->nb_col)
	    {
	      zmatbeg[i]= zmatbeg[i-1] + 1;
	      for(j=0;j<qqp->n;j++)
		{
		  if (!Zero_BB(qqp->nb_var_y[ij2k(l,j,qqp->n)]))
		    zmatbeg[i]= zmatbeg[i] + 2;
		}
	      i++;
	    }
	}
    }
    
  int numnz;
  int cpt,cpt_bis;
  /* t variable appears in 1 constraint if none y_ij is created for a fixed i, and 2 more constraints for each j where y_ij is created*/
 
  numnz=zmatbeg[qqp->nb_col-1] + 1;
  for(j=0;j<qqp->n;j++)
    if (!Zero_BB(qqp->nb_var_y[ij2k(qqp->nb_int-1,j,qqp->n)]))
      numnz=numnz + 2;
    
  for(i=0;i<qqp->nb_col-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  zmatcnt[qqp->nb_col - 1] = numnz - zmatbeg[qqp->nb_col-1];
  
 
  /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
  /* Contraints over x (order x_1 x_2 ... x_n ) */
  /* Step 1 : from x_1 to x_nb_int*/
  for(i=0;i<qqp->ind_max_int;i++)
    {
      /*Ax = b*/
      for(j= 0;j<qqp->ind_eg;j++)
	{ 
	  zmatval[s]=qqp->a[ij2k(j,i,qqp->n)];
	  zmatind[s]=qqp->cont[j];
	  s++;
	}
       
      /* Dx =e*/
      for(j=qqp->ind_eg;j<qqp->ind_ineg;j++)
	{
	  zmatval[s]=qqp->d[ij2k(j-qqp->ind_eg,i,qqp->n)];
	  zmatind[s]=qqp->cont[j];
	  s++;
	}
        /*Aq,y + Aq_0 x = bq*/
      for(j=qqp->ind_ineg;j<qqp->ind_eg_q;j++)
	{ 
	  zmatval[s]=2*qqp->aq[ijk2l(j-qqp->ind_ineg,0,i+1,qqp->n+1,qqp->n+1)];
	  zmatind[s]=j;
	  s++;
	  
	}
       
      /*Dq,y + Dq_0 x <= eq*/
      for(j=qqp->ind_eg_q;j<qqp->ind_ineg_q;j++)
	{
	  zmatval[s]=2*qqp->dq[ijk2l(j-qqp->ind_eg_q,0,i+1,qqp->n+1,qqp->n+1)];
	  zmatind[s]=j;
	  s++;
	} 
	
          
      /* x_i = sum_k 2^k t_ik*/
      zmatind[s]=qqp->ind_ineg_q+i;
      zmatval[s]=1;  
      s++;  
      
        
      /* z_ijk <= x_j*/
      
      cpt=0;
      cpt_bis=0;
      for (j=0;j<qqp->nb_int;j++)
	{
	  lim=log_base_2(qqp->u[j])+1;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		zmatind[s]=qqp->cont[qqp->ind_z_1 + i*lim + cpt_bis + cpt + k];
		zmatval[s]=-1;
		s++;
	      }
	  cpt=cpt+ (qqp->n - qqp->nb_int)*lim;
	  cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
	}
      
      /* z_ijk >= x_j + u_jt_ik - u_j*/
     
      cpt=0;
      cpt_bis=0;
      for (j=0;j<qqp->nb_int;j++)
	{
	  lim=log_base_2(qqp->u[j])+1;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		zmatind[s]=qqp->cont[qqp->ind_z_2 + i*lim + cpt_bis+ cpt+k];
		zmatval[s]=1;
		s++;
	    }
	  cpt=cpt+(qqp->n- qqp->nb_int)*lim;
	  cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
	}
      
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<qqp->ind_max_x;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    {
	      if(j>i)
		{ 
		  zmatind[s] =qqp->cont[qqp->ind_z_4+ i* (2*(qqp->n) -i +1)/2 + j-i];
		  zmatval[s]=qqp->u[j]; 
		  s++;
		}
	      if(i>j)
		{ 
		  zmatind[s]=qqp->cont[qqp->ind_z_4 + j* (2*(qqp->n) -j +1)/2  + i-j];
		  zmatval[s]=qqp->u[j];
		  s++;
		}   
	      if(i==j)
		{
		  zmatind[s]=qqp->cont[qqp->ind_z_4  + j* (2*(qqp->n) -j +1)/2];
		  zmatval[s]=2*qqp->u[i];
		  s++;
		}
	    }
	}
       
      /* y_ii >= x_i*/
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,i,qqp->n)]))
	{
	  zmatind[s]=qqp->cont[qqp->ind_y_1 + i];
	  zmatval[s]=1;
	  s++;
	}
    }
  
        
  /* Step 2 : from x_nb_int+1 to x_n*/
  for(i=qqp->ind_max_int;i<qqp->ind_max_x;i++)
    {
      /*Ax = b*/
      if (i<qqp->ind_max_x)
	for(j= 0;j<qqp->ind_eg;j++)
	  { 
	    zmatval[s]=qqp->a[ij2k(j,i,qqp->n)];
	    zmatind[s]=qqp->cont[j];
	    s++;
	  }
      
      /* Dx +s =e*/
      for(j=qqp->ind_eg;j<qqp->ind_ineg;j++)
	{
	  zmatval[s]=qqp->d[ij2k(j-qqp->ind_eg,i,qqp->n)];
	  zmatind[s]=qqp->cont[j];
	  s++;
	}
       
      /*Aq,y + Aq_0 x = bq*/
      if (i<qqp->ind_max_x)
	for(j=qqp->ind_ineg;j<qqp->ind_eg_q;j++)
	  { 
	    zmatval[s]=2*qqp->aq[ijk2l(j-qqp->ind_ineg,0,i+1,qqp->n+1,qqp->n+1)];
	    zmatind[s]=j;
	    s++;
	  }
      
      /*Dq,y + Dq_0 x  + sq = eq*/
      for(j=qqp->ind_eg_q;j<qqp->ind_ineg_q;j++)
	  {
	    if (i<qqp->ind_max_x)
	      {
		zmatval[s]=2*qqp->dq[ijk2l(j-qqp->ind_eg_q,0,i+1,qqp->n+1,qqp->n+1)];
		zmatind[s]=j;
	      }
	    else
	      {
		zmatval[s]=0;
		zmatind[s]=j;
	      }
	    s++;
	  } 
      
      /* z_jik <= x_i*/
      cpt =0;
      cpt_bis=0;
      for (j=0;j<qqp->nb_int;j++)
	{
	  lim=log_base_2(qqp->u[j])+1;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		zmatind[s]=qqp->cont[qqp->ind_z_1 + i*lim + cpt_bis + cpt + k];
		zmatval[s]=-1;
		s++;
	      }
	  cpt=cpt+(qqp->n - qqp->nb_int)*lim;
	  cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
	}
      
      /* z_jik >= x_i + u_it_jk - u_i*/
      cpt =0;
      cpt_bis=0;
      for (j=0;j<qqp->nb_int;j++)
	{
	  lim=log_base_2(qqp->u[j])+1;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		zmatind[s]=qqp->cont[qqp->ind_z_2 + i*lim + cpt_bis+cpt+ k];
		zmatval[s]=1;
		s++;
	      }
	  cpt=cpt+(qqp->n - qqp->nb_int)*lim;
	  cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
	}
       
       
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<qqp->ind_max_int;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	  {	  
	    zmatind[s]=qqp->cont[qqp->ind_z_4 + j* (2*(qqp->n ) -j +1)/2  + i-j];
	    zmatval[s]=qqp->u[j];
	    s++;
	  }
    }   

  /* Constraints on y in order y_11 y_12 .. y_1n , ... ,y_n1 .. y_nn */   
  /* Step 1 : from y_11 to y_nb_int n*/
  for(i=0;i<qqp->ind_max_int;i++)
    for(j=0;j<qqp->ind_max_x;j++)
      {
	/*Aq,y + Aq_0 x = bq*/
	for(k=qqp->ind_ineg;k<qqp->ind_eg_q;k++)
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	     { 
	       if (i==j)
		 zmatval[s]=qqp->aq[ijk2l(k-qqp->ind_ineg,i+1,j+1,qqp->n+1,qqp->n+1)];
	       else
		 {
		   if(i<j)
		     zmatval[s]=2*qqp->aq[ijk2l(k-qqp->ind_ineg,i+1,j+1,qqp->n+1,qqp->n+1)];
		 }
	       zmatind[s]=k;
	       s++;
	     }

	/* Dq,y + Dq_0 x <=eq*/
	for(k=qqp->ind_eg_q;k<qqp->ind_ineg_q;k++)
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	     {
	       if (i==j)
		 zmatval[s]=qqp->dq[ijk2l(k-qqp->ind_eg_q,i+1,j+1,qqp->n+1,qqp->n+1)];
	       else
		 {
		   if(i<j)
		     zmatval[s]=2*qqp->dq[ijk2l(k-qqp->ind_eg_q,i+1,j+1,qqp->n+1,qqp->n+1)];
		 }
	       zmatind[s]=k;
	       s++;
	     } 
	
	/* y_ij = sum_k 2^k z_ijk*/
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	  {
	    zmatind[s] =qqp->cont[qqp->ind_dec_x+ i*(qqp->n ) +j];
	    zmatval[s]=1; 
	    s++;
	  }
	
	if(i==j)
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    {
	      
	      /*y_ii >= 2u_ix_i - u_i^2*/ 
	      zmatind[s] =qqp->cont[qqp->ind_z_4 + i* (2*(qqp->n  ) -i +1)/2];
	      zmatval[s]=-1; 
	      s++;   
	      
	      /*y_ii >= x_i*/ 
	      zmatind[s] =qqp->cont[qqp->ind_y_1+i];
	      zmatval[s]=-1; 
	      s++;
	    }
	 
	if(j>i)
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    {	  
	      /*y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
	      zmatind[s] =qqp->cont[qqp->ind_z_4+ i* (2*(qqp->n  ) -i +1)/2 + j-i];
	      zmatval[s]=-1; 
	      s++;
	      
	      /*y_ij = y_ji*/ 
	      zmatind[s] =qqp->cont[qqp->ind_y_2+ i* (2*(qqp->n  ) -i-1)/2 + j-i-1];
	      zmatval[s]=1; 
	      s++;
	    }
	 
	if(i>j)
	  /*y_ij = y_ji*/   
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    {
	      zmatind[s] =qqp->cont[qqp->ind_y_2+ j* (2*(qqp->n ) -j-1)/2 + i-j-1];
	      zmatval[s]=-1; 
	      s++;
	    }
      }
  
  /* Step 2 : from y_nb_int 1 to y_n nb_int*/ 
  for(i=qqp->ind_max_int;i<qqp->ind_max_x;i++)
    for(j=0;j<qqp->ind_max_int;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	  /*y_ij = y_ji*/   
	  zmatind[s] =qqp->cont[qqp->ind_y_2+ j* (2*(qqp->n) -j-1)/2 + i-j-1];
	  zmatval[s]=-1; 
	  s++;
	}
   
  /* constraints over z  (order z_110 z_111 .. z_nn0 z_nnk) */    
  cpt=0;
  cpt_bis=0;
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->n;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		/* y_il = sum_k 2^k z_ijk */
		zmatind[s] =qqp->cont[qqp->ind_dec_x+ j + i*(qqp->n)];         
		zmatval[s]= -pow(2,k);
		s++;  
		
		/* z_ijk <= u_j_t_ik */
		zmatind[s] = qqp->cont[qqp->ind_dec_y + cpt_bis + j*lim + k+cpt];
		zmatval[s]=1; 
		s++;
		
		/* z_ijk <= x_j */
		zmatind[s] =qqp->cont[ qqp->ind_z_1+ cpt_bis+ j*lim + k+cpt];
		zmatval[s]=1; 
		s++;
		
		/* z_ijk >= x_j + u_j_t_ik - u_j */
		zmatind[s] = qqp->cont[qqp->ind_z_2 + cpt_bis+j*lim + k+cpt];
		zmatval[s]=-1; 
		s++;
		
		/* z_ijk >= 0 */
		zmatind[s] = qqp->cont[qqp->ind_z_3 + cpt_bis+j*lim + k+cpt];
		zmatval[s]=-1; 
		s++;
	      }
	}
      cpt=cpt+(qqp->n  - qqp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
    }

   
  /* constraints over t (order t_10 t_11 .. t_n0 t_nk) */
  cpt=0;
  cpt_bis=0;
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<lim;j++)
	{
	  /* x_i = sum_k 2^k t_ik */
	  zmatind[s] = qqp->ind_ineg_q + i;         
	  zmatval[s]= -pow(2,j);
	  s++;  
	   
	  for(k=0;k<qqp->n;k++)
	    if (!Zero_BB(qqp->nb_var_y[ij2k(i,k,qqp->n)]))
	      {
		/* z_ikj <= u_k t_ij */
		zmatind[s] = qqp->cont[qqp->ind_dec_y + j+ cpt_bis+ k*lim+cpt];
		zmatval[s]=-qqp->u[k]; 
		s++;
	      }
	   
	  for(k=0;k<qqp->n;k++)
	    if (!Zero_BB(qqp->nb_var_y[ij2k(i,k,qqp->n)]))
	      {
		/* z_ikj >= x_k + u_k t_ij - u_k */
		zmatind[s] = qqp->cont[qqp->ind_z_2 + j + cpt_bis + k*lim+cpt];
		zmatval[s]=qqp->u[k]; 
		s++;
	      }
	}
      cpt=cpt+(qqp->n - qqp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
    }
  
 
  /* Second member of constrains */ 
  zrhs = alloc_vector_d (qqp->nb_row);
  i=0;
   
  /* Ax = b */
  for(j=0;j<qqp->m;j++)
    {
      zrhs[i]=qqp->b[j];
      i++;
    }
   
  /* Dx <= e*/
  for(j=0;j<qqp->p;j++)
    {
      zrhs[i]=qqp->e[j];
      i++;
    }
    
  /*Aq,y + Aq_0 x = bq*/
  for(j=0;j<qqp->mq;j++)
    {
      zrhs[i]=qqp->bq[j];
      i++;
    }
  
  /*Dq,y + Dq_0 x <= eq*/  
  for(j=0;j<qqp->pq;j++)
    {
      zrhs[i]=qqp->eq[j];
      i++;
    }
  for(;i<qqp->ind_dec_x+ qqp->nb_dec_y + qqp->nb_z_1 + qqp->nb_z_2;i++)
    {
      zrhs[i]=0;
    }
  
  int ind=i;
  for(i=0;i<qqp->nb_int ;i++)
    {
      for(j=0;j<qqp->n;j++)
	{
	  lim=log_base_2(qqp->u[i])+1;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	    {
	      zrhs[ind]=qqp->u[j];
	      ind++;
	    }
	}
    }
   
  i=ind;
  for(;i<qqp->ind_dec_x+ qqp->nb_dec_y + qqp->nb_z_1 + qqp->nb_z_2 + qqp->nb_z_3  + qqp->nb_z_4 ;i++)
    {
      zrhs[i]=0;
    }
   
  for(j=0;j<qqp->nb_int;j++)
    {
      for(k=j;k<qqp->n ;k++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)]))
	  {
	    zrhs[i]= qqp->u[j]*qqp->u[k];
	    i++;
	  }
    }
   
  for(;i<qqp->nb_row;i++)
    zrhs[i]=0;
  
  /* Sense of constraints */
  zsense = alloc_string(qqp->nb_row);

  for(i=0;i<qqp->ind_eg;i++)
    zsense[i]= 'E';
  for(;i<qqp->ind_ineg;i++)
    zsense[i]= 'L';
  for(i;i<qqp->ind_eg_q;i++)
    zsense[i]= 'E';
  for(;i<qqp->ind_ineg_q;i++)
    zsense[i]= 'L';
  for(;i<qqp->ind_dec_x+ qqp->nb_dec_y;i++)
    zsense[i]= 'E';
  for(;i<qqp->ind_dec_x+ qqp->nb_dec_y + qqp->nb_z_1 + qqp->nb_z_2+ qqp->nb_z_3+ qqp->nb_z_4+ qqp->nb_y_1+ qqp->nb_y_2;i++)
    zsense[i]= 'L';
  for(;i<qqp->nb_row;i++)
    zsense[i]= 'E';
    

  
 
  /***********************************************************************************/
  /****************************** Bounds on variables ********************************/
  /***********************************************************************************/
   
  zlb= alloc_vector_d(qqp->nb_col);
  for(i=0;i<qqp->n;i++)
    zlb[i]=qqp->l[i];

  k=qqp->n;
  for(i=0;i<qqp->nb_int;i++)
    for(j=0;j<qqp->n ;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	  zlb[k]=qqp->l[i]*qqp->l[j];
	  k++;
	}
    
  for(i=qqp->nb_int;i<qqp->n ;i++)
    for(j=0;j<qqp->nb_int;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	  zlb[k]=qqp->l[i]*qqp->l[j];
	  k++;
	}
  for(;k<qqp->ind_max_x+qqp->nb_y+qqp->nb_z;k++)
    zlb[k]=-CPX_INFBOUND;
  for(;k<qqp->nb_col;k++)
    zlb[k]=0;
  
  zub = alloc_vector_d(qqp->nb_col);
   
  for(i=0;i<qqp->n ;i++)
    zub[i]=qqp->u[i];
   
  k=qqp->n;
  for(i=0;i<qqp->nb_int;i++)
    for(j=0;j<qqp->n ;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	  zub[k]=qqp->u[i]*qqp->u[j];
	  k++;
	}
    
  for(i=qqp->nb_int;i<qqp->n ;i++)
    for(j=0;j<qqp->nb_int;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	  zub[k]=qqp->u[i]*qqp->u[j];
	  k++;
	}
   
  for(;k<qqp->ind_max_x+qqp->nb_y+qqp->nb_z;k++)
    zub[k]=CPX_INFBOUND;
  for(;k<qqp->nb_col;k++)
    zub[k]=1;

  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
 
  zobj = alloc_vector_d (qqp->nb_col);
  for(i=0;i<qqp->ind_max_x;i++)
    {
      if(i<qqp->n)
	    zobj[i]= qqp->c[i];
      else
	zobj[i]=0;
    }



  for(i=qqp->ind_max_x;i<qqp->nb_col;i++)
    zobj[i]=0;
   
  i=qqp->ind_max_x;
  for(j=0;j<qqp->ind_max_int;j++)
    for(k=0;k<qqp->n ;k++)
      {
	if ((j==k) && (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)])))
	  {
	    zobj[i]= qqp->q[ij2k(j,k,qqp->n)] - qqp->new_q[ij2k(j,k,qqp->n)]/2 + qqp->new_lambda_min; 
	    i++;
	  }
	else
	  if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)]))
	    {
	      zobj[i]= qqp->q[ij2k(j,k,qqp->n)] - (qqp->new_q[ij2k(j,k,qqp->n)])/2;
	      i++;
	    }
      }
   
  for(j=qqp->ind_max_int;j<qqp->ind_max_x;j++)
    for(k=0;k<qqp->ind_max_int;k++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)]))
	{
	  zobj[i]= qqp->q[ij2k(j,k,qqp->n)] - (qqp->new_q[ij2k(j,k,qqp->n)])/2; 
	  i++;
	}

  int numqnz = nb_non_zero_matrix(qqp->new_q, qqp->n,qqp->n);

  zqmatbeg  = alloc_vector(qqp->nb_col);
  zqmatbeg[0]=0;
  for (i=1;i<=qqp->ind_max_x;i++)
    zqmatbeg[i]=zqmatbeg[i-1]+ qqp->ind_max_x;
   
  for (i;i<qqp->nb_col;i++)
    zqmatbeg[i]=zqmatbeg[qqp->ind_max_x];
   
  zqmatcnt   = alloc_vector(qqp->nb_col);
  for (i=0;i<qqp->ind_max_x;i++)
    zqmatcnt[i] = zqmatbeg[i+1] - zqmatbeg[i];
   
  for (i;i<qqp->nb_col;i++)
    zqmatcnt[i]=0; 

  zqmatind   = alloc_vector(qqp->ind_max_x*qqp->ind_max_x);
  zqmatval   = alloc_vector_d(qqp->ind_max_x*qqp->ind_max_x);
   
  k=0;
  for(i=0;i<qqp->ind_max_x;i++)
    for(j=0;j<qqp->ind_max_x;j++)
      { 
	if (i==j)
	  zqmatval[k] = qqp->new_q[ij2k(i,j,qqp->n)] - ((qqp->new_lambda_min)*2); 
	else
	  zqmatval[k] =  qqp->new_q[ij2k(i,j,qqp->n)]; 
	 
	zqmatind[k] = j;
	k++;
      }
 
  zctype = alloc_string(qqp->nb_col); 
   
    
  for(i=0;i<qqp->ind_max_x;i++)
    zctype[i]= 'C';
  
  for(k=0;k<qqp->nb_int;k++)
    for(l=0;l<qqp->ind_max_x;l++)
   	if (!Zero_BB(qqp->nb_var_y[ij2k(k,l,qqp->n)]))
	  {
	    zctype[i]='C';
	    i++;
	  }
   
  for(k=qqp->nb_int;k<qqp->ind_max_x;k++)
    for(l=0;l<qqp->nb_int;l++)
      	if (!Zero_BB(qqp->nb_var_y[ij2k(k,l,qqp->n)]))
	  {
	    zctype[i]='C';
	    i++;
	  }
  
  cpt=0;
  cpt_bis=0;
  for(j=0;j<qqp->nb_int;j++)
    {
      lim=log_base_2(qqp->u[j])+1;
      for(l=0;l<qqp->n;l++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(j,l,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		 zctype[i]= 'C';
		 i++;
	      }
	}
      cpt=cpt+(qqp->n  - qqp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qqp->nb_int)*lim;
    }

  for(i;i<qqp->nb_col;i++)
    zctype[i]= 'B';
  
   strcpy (zprobname, "miqcr_cplex");
   
  *numcols_p   = qqp->nb_col;
  *numrows_p   = qqp->nb_row;
  *objsen_p    = CPX_MIN;   /* The problem is minimization */
   
  *probname_p  = zprobname;
  *obj_p       = zobj;
  *rhs_p       = zrhs;
  *sense_p     = zsense;
  *matbeg_p    = zmatbeg;
  *matcnt_p    = zmatcnt;
  *matind_p    = zmatind;
  *matval_p    = zmatval;
  *qmatbeg_p    = zqmatbeg;
  *qmatcnt_p    = zqmatcnt;
  *qmatind_p    = zqmatind;
  *qmatval_p    = zqmatval;
  *lb_p        = zlb;
  *ub_p        = zub;
  *ctype_p = zctype; 
  

  
   
  return (status);
} /* END setproblemdata */

void miqcr_cplex()
{
  char     *probname = NULL;  
  int      numcols;
  int      numrows;
  int      objsen;
  double   *obj = NULL;
  double   *rhs = NULL;
  char     *sense = NULL;
  int      *matbeg = NULL;
  int      *matcnt = NULL;
  int      *matind = NULL;
  double   *matval = NULL;
  double   *lb = NULL;
  double   *ub = NULL;
  char     *ctype = NULL;
  int      *qmatbeg = NULL;
  int      *qmatcnt = NULL;
  int      *qmatind = NULL;
  double   *qmatval = NULL;
  
 

  int      solstat,i,j;
  
  qqp->sol = alloc_vector_d (qqp->nb_col);

  /*+1 for constant term*/

  CPXENVptr    env = NULL;
  CPXLPptr      lp = NULL;
  int           status;
   
  int           cur_numrows, cur_numcols;
  int file_sol;
  int stdoutsave;
  
 
  
  /* Initialize the CPLEX environment */
  int k,l,temp;
    
  /* load program into cplex*/
  env = CPXopenCPLEX (&status);
  char  errmsg[1024];
  if ( env == NULL ) 
    {
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return;
    }
   
  if ( env == NULL ) {
    printf ("Could not open CPLEX environment.\n");
    return;
  }
   
  /* Turn on output to the screen */

  status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  if ( status ) 
    {
      printf ("Failure to turn on screen indicator, error %d.\n", status);
      return;
    }
   
  /* use quadratic solver (cplex12.6)*/

  status = CPXsetintparam (env, 4012, 0);
  if ( status )
    {
      printf ("Failure to turn on screen indicator, error %d.\n", status);
      return;
    }

  /*var selection strategie*/
  status = CPXsetintparam (env, CPX_PARAM_VARSEL, VARSEL);
  if ( status )
    {
      printf ("Failure to change var selection strategie, error %d.\n", status);
      return;
    }

  status = CPXsetdblparam (env, CPX_PARAM_EPGAP,REL_GAP);
  if ( status )
    {
      printf ("Failure to change relative gap, error %d.\n", status);
      return;
    }

  if (qqp->n == qqp->nb_int)
    {
      status = CPXsetdblparam (env, CPX_PARAM_EPAGAP,ABS_GAP);
      if ( status ) 
	{
	  printf ("Failure to change absolute gap, error %d.\n", status);
	  return;
	}
    }
  
  status = CPXsetdblparam (env, CPX_PARAM_TILIM,(double)TIME_LIMIT_CPLEX);
  
  if ( status ) 
    {
      printf ("Failure to change time limit, error %d.\n", status);
      return;
    }

  /* /\* Set memory parameter*\/ */
  status = CPXsetdblparam (env, CPX_PARAM_WORKMEM,MEMORY);
  if ( status )
    {
      printf ("Failure to change memory, error %d.\n", status);
      return;
    }

   
  /* Fill in the data for the problem.  */

  status = setqpproblemdata_miqcr_cplex ( &probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval,&qmatbeg, &qmatcnt, &qmatind, &qmatval, &lb, &ub, &ctype);
   
  if ( status ) 
    {
      printf ( "Failed to build problem data arrays.\n");
      return;
    }

  /* Create the problem. */
  lp = CPXcreateprob (env, &status, probname);

  if ( lp == NULL ) 
    {
      printf ("Failed to create LP.\n");
      return;
    }
   
  /* Now copy the problem data into the lp */
  status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs, 
		      sense, matbeg, matcnt, matind, matval,
		      lb, ub, NULL);
   
  if ( status )
    {
      printf ( "Failed to copy problem data.\n");
      return ;
    }

  /* Now copy the ctype array */
  status = CPXcopyctype (env, lp, ctype);
  if ( status )
    {
      printf ( "Failed to copy ctype\n");
      return;
    }

  /* Now add quadratic part of the objective function */
  status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
  if ( status ) 
    {
      printf ("Failed to copy quadratic matrix.\n");
      return;
    }
 
  /* Optimize the problem and obtain solution. */
  status = CPXmipopt (env, lp);
  if ( status ) 
    {
      printf ( "Failed to optimize MIQP.\n");
      return;
    }
    
  solstat = CPXgetstat (env, lp);
  
  /* Write the output to the screen. */
  status = CPXgetmipobjval (env, lp, &ub_bb);
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return;
    }
   
  /*Print the Optimal solution into the save file */  
     
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  

  
  status = CPXgetmipx (env, lp, best_x, 0, cur_numcols-1);
  if ( status ) 
    {
      printf ( "Failed to get optimal integer x.\n");
      return;
    }
       
  /*number of branch-and-bound nodes*/
  nodecount = CPXgetnodecnt (env, lp); 
  ub_bb = ub_bb + qqp->cons;
 
  
 

  /* Free up the problem as allocated by CPXcreateprob, if necessary */
  if ( lp != NULL )  
    {
      status = CPXfreeprob (env, &lp);
      if ( status )
	{
	  printf ( "CPXfreeprob failed, error code %d.\n", status);
	  return;
	}
    }
   
  /* Free up the CPLEX environment, if necessary */

  if ( env != NULL ) 
    {
      status = CPXcloseCPLEX (&env);
        
      if ( status ) 
	{
	  fprintf (stderr, "Could not close CPLEX environment.\n");
	  CPXgeterrorstring (env, status, errmsg);
	  fprintf (stderr, "%s", errmsg);
	  return;
	}
    }
   
  /* Free up the problem data arrays, if necessary. */

  free_and_null ((char **) &probname);
  free_and_null ((char **) &obj);
  free_and_null ((char **) &rhs);
  free_and_null ((char **) &sense);
  free_and_null ((char **) &matbeg);
  free_and_null ((char **) &matcnt);
  free_and_null ((char **) &matind);
  free_and_null ((char **) &matval);
  free_and_null ((char **) &qmatbeg);
  free_and_null ((char **) &qmatcnt);
  free_and_null ((char **) &qmatind);
  free_and_null ((char **) &qmatval);
  free_and_null ((char **) &lb);
  free_and_null ((char **) &ub);
  free_and_null ((char **) &ctype);

  
   
}
