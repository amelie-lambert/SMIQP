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
#include <ctype.h>

#include <ilcplex/cplex.h>

#include"obbt.h"
#include "utilities.h"

extern MIQCP qp;
extern C_MIQCP cqp;

/* parameters for CPLEX */
extern int MEMORY;

extern double EPS_OBBT;
extern int OBBT_QUAD;
extern int OBBT_CONT;

int setquadconst_cplex_obbt ( double upperbound,int **linnzcnt, int **quadnzcnt, double **qrhs, char **qsense, int **linind, double **linval,int **qconrow_p, int **qconcol_p, double **qconval_p)
{
  int status =0;
  int      *zlinnzcnt= NULL;
  int      *zquadnzcnt= NULL;
  double   *zqrhs= NULL;
  char     *zqsense= NULL;
  int      *zlinind= NULL;
  double   *zlinval= NULL;
  int      *zqconrow = NULL;
  int      *zqconcol = NULL;
  double   *zqconval = NULL;

  int i,j,s;

  
  /***********************************************************************************/
  /****************************** Quadratic Constraints ******************************/
  /*********  x^t S^* x + <Q- S^* ,y> + c^T x <= upper_bound
  /***********************************************************************************/
  int ind_lin=0;
  int ind_q1=0;
  int ind_q2=0;
  int r;
  int nb_quad_cont = 1; 

  /* non zero elements + second member*/
  zlinnzcnt = alloc_vector(nb_quad_cont);
  zquadnzcnt= alloc_vector(nb_quad_cont);
  zqrhs= alloc_vector_d(nb_quad_cont);
  zqsense=alloc_string(nb_quad_cont);
  // + 1 for the constant term 
  zlinnzcnt[0]=cqp->n + (cqp->n +1)*cqp->n/2 + 1;
  zquadnzcnt[0]=cqp->n*cqp->n;
  zqrhs[0]=upperbound + EPS_OBBT;
  zqsense[0]='L';
    
    
  
    
  /* linear terms*/
  int num_lin_terms =zlinnzcnt[0];
  
  zlinind = alloc_vector(num_lin_terms);
  zlinval = alloc_vector_d(num_lin_terms);

   
  /*x^TAq x + Aq_0 x = bq*/
  s=0;
  /* x variables */
  for (i=0;i<cqp->n;i++)
    {
      zlinval[s] = cqp->c[i];
      zlinind[s] = i;
      s++;
    }
  /* y variables */
  for (i=0;i<cqp->n;i++)
    for (j=i;j<cqp->n;j++)
      {
	if (i==j)
	  {
	    zlinval[s] = cqp->q[ij2k(i,j,cqp->n)] - cqp->new_q[ij2k(i,j,cqp->n)]/2 + cqp->new_lambda_min;
	    zlinind[s] = cqp->n + i* (2*(cqp->n) -i +1)/2 + j-i;
	    s++;
	  }
	else
	  {
	    zlinval[s] =2*(cqp->q[ij2k(i,j,cqp->n)] - (cqp->new_q[ij2k(i,j,cqp->n)])/2);
	    zlinind[s] =cqp->n + i* (2*(cqp->n) -i +1)/2 + j-i;
	    s++;
	  }
      }
  // constant term 
  zlinval[s] = 1;
  zlinind[s] = cqp->n +  (cqp->n +1)*cqp->n/2 ;
  
   
  /* quadratic terms */
  int num_quad_terms=zquadnzcnt[0];

  
  zqconrow = alloc_vector(num_quad_terms);
  zqconcol = alloc_vector(num_quad_terms);
  zqconval= alloc_vector_d(num_quad_terms);
   
  /*x^TAq x + Aq_0 x = bq*/
  s=0;
  for (i=0;i<cqp->n;i++)
    for (j=0;j<cqp->n;j++)
      {
	if (i==j)
	  zqconval[s] = cqp->new_q[ij2k(i,j,cqp->n)]/2 - (cqp->new_lambda_min) ;
	else
	  zqconval[s] = cqp->new_q[ij2k(i,j,cqp->n)]/2;
	zqconrow[s] = i; 
	zqconcol[s] = j; 
	s++;
      }


  
  
  *linnzcnt = zlinnzcnt;
  *quadnzcnt = zquadnzcnt;
  *qrhs = zqrhs;
  *qsense = zqsense;
  *linind = zlinind;
  *linval = zlinval;
  *qconrow_p   = zqconrow;
  *qconcol_p   = zqconcol;
  *qconval_p   = zqconval;
 
   return (status);
}

int setcqpproblemdata_cplex_obbt_var ( char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, double **lb_p, double **ub_p, char ** ctype_p, int ind_var,int sense_opt)
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
   
  int i,j,k,l,s;

  double * new_q;

  zprobname = alloc_string(16); 
  
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
  // + 1 variable for the constant term
  int nbcols = cqp->n + (cqp->n)*(cqp->n+1)/2 + 1;
  zmatbeg   = alloc_vector(nbcols);   
  zmatcnt   = alloc_vector(nbcols); 

 
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i<=cqp->n;i++)
    zmatbeg[i]= zmatbeg[i-1] + 4*(cqp->n-1) + 3 + cqp->m+cqp->p +cqp->mq +cqp->pq;
    
    
  /* y  */
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
      {
	zmatbeg[i]=zmatbeg[i-1];
	if (j==k)
	  zmatbeg[i]= zmatbeg[i] +3;
	else
	  zmatbeg[i]= zmatbeg[i] + 4;
	
	if (!Zero_BB(cqp->mq))
	  zmatbeg[i]= zmatbeg[i] + cqp->mq;
	if (!Zero_BB(cqp->pq))
	  zmatbeg[i]= zmatbeg[i] + cqp->pq;
	i++;
      }
      
 
 
  int numnz;
  
  numnz = zmatbeg[nbcols-1];
  
  for(i=0;i<nbcols-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  zmatcnt[nbcols - 1] = 0;
  

  
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
  
  /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  int ind_eg = cqp->m;
  int ind_ineg = ind_eg+cqp->p;
  int ind_eg_q = ind_ineg + cqp->mq;
  int ind_ineg_q = ind_eg_q + cqp->pq;
  int ind_y_1  = ind_ineg_q + (cqp->n+1)*cqp->n/2;
  int ind_y_2 = ind_y_1 + (cqp->n-1)*cqp->n/2 ;
  int ind_y_3 = ind_y_2 + (cqp->n+1)*cqp->n/2 ;
  int nbrows = ind_y_3 + (cqp->n+1)*cqp->n/2 ;
 
  /* Contraints on x in order x_1 x_2 ... x_n */
  for(i=0;i<cqp->n;i++)
    {
      
      /*Ax = b*/
      for(j= 0;j<ind_eg;j++)
	{ 
	  zmatval[s]=cqp->a[ij2k(j,i,cqp->n)];
	  zmatind[s]=j;
	  s++;
	}
       
      /* Dx =e*/
      for(j=ind_eg;j<ind_ineg;j++)
	{
	  zmatval[s]=cqp->d[ij2k(j-ind_eg,i,cqp->n)];
	  zmatind[s]=j;
	  s++;
	}
      
      /*Aq,y + Aq_0 x = bq*/
      for(j=ind_ineg;j<ind_eg_q;j++)
      	{
      	  zmatval[s]=2*cqp->aq[ijk2l(j-ind_ineg,0,i+1,cqp->n+1,cqp->n+1)];
      	  zmatind[s]=j;
      	  s++;
      	}
      
      /*Dq,y + Dq_0 x <= eq*/
      for(j=ind_eg_q;j<ind_ineg_q;j++)
	{
	  zmatval[s]=2*cqp->dq[ijk2l(j- ind_eg_q,0,i+1,cqp->n+1,cqp->n+1)];
	  zmatind[s]=j;
	  s++;
	} 
        
      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
      for (j=0;j<cqp->n;j++)
	{
	  if(j>i)
	      { 
		zmatind[s]= ind_ineg_q + i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=-(double)cqp->u[j];
		s++;
	      }
	  if(i>j)
	    { 
		zmatind[s]=ind_ineg_q + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=-(double)cqp->l[j];
		s++;
	    } 
	  if(i==j)
	    {
		zmatind[s]=ind_ineg_q + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=-((double)cqp->l[i]+(double)cqp->u[j]);
		s++;
	      }
	}
      
	        
      /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
      for (j=0;j<cqp->n;j++)
	{
	  if(j>i)
	    { 
	      zmatind[s]= ind_y_1 + i* (2*(cqp->n) - i -1)/2 + j-i-1;
	      zmatval[s]=-(double)cqp->l[j];
	      s++;
	    }
	  if(i>j)
	    { 
	      zmatind[s]=ind_y_1 + j* (2*(cqp->n) -j+ -1)/2  + i-j-1;
	      zmatval[s]=-(double)cqp->u[j];
	      s++;
	    } 
	  }
      
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<cqp->n;j++)
	{
	    if(j>i)
	      { 
		zmatind[s] =ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)cqp->u[j]; 
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=ind_y_2 + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=(double)cqp->u[j];
	      s++;
	      }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_2  + j* (2*(cqp->n) -j +1)/2;
		zmatval[s]=2*(double)cqp->u[i];
		s++;
	      }
	  }
      
      /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
      for (j=0;j<cqp->n;j++)
	  {
	    if(j>i)
	      { 
		zmatind[s] =ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)cqp->l[j]; 
		s++;
	      }
	    if(i>j)
	    { 
	      zmatind[s]=ind_y_3 + j* (2*(cqp->n) -j +1)/2  + i-j;
	      zmatval[s]=(double)cqp->l[j];
	      s++;
	    }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_3  + j* (2*(cqp->n) -j +1)/2;
		zmatval[s]=2*(double)cqp->l[i];
		s++;
	      }
	  }
    }

      
  /* Constraints on y in order y_11 y_12 .. y_1n , ... ,y_n1 .. y_nn */   
  
  for(i=0;i<cqp->n;i++)
    for(j=i;j<cqp->n;j++)
      {
	  
	/*Aq,y + Aq_0 x = bq*/
	for(k=ind_ineg;k<ind_eg_q;k++)
	  {
	    if (i==j)
	      {
		zmatval[s]=cqp->aq[ijk2l(k-ind_ineg,i+1,j+1,cqp->n+1,cqp->n+1)];
		zmatind[s]=k;
		s++;
	  	}
	    else
	      {
	  	  zmatval[s]=2*cqp->aq[ijk2l(k-ind_ineg,i+1,j+1,cqp->n+1,cqp->n+1)];
	  	  zmatind[s]=k;
	  	  s++;
	      }
	  }
	  
	  /* Dq,y + Dq_0 x <=eq*/
	  for(k=ind_eg_q;k<ind_ineg_q;k++)
	    {
	      if (i==j)
		{
		  zmatval[s]=cqp->dq[ijk2l(k-ind_eg_q,i+1,j+1,cqp->n+1,cqp->n+1)];
		  zmatind[s]=k;
		  s++;
		}
	      else
		{
		  zmatval[s]=2*cqp->dq[ijk2l(k-ind_eg_q,i+1,j+1,cqp->n+1,cqp->n+1)];
		  zmatind[s]=k;
		  s++;
		}
	    }

 
	  if(i==j)
	    {	  
	      /*y_ii <= l_ix_i+u_ix_i - l_iu_i*/ 
	      zmatind[s] =ind_ineg_q + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ii >= 2u_ix_i - u_i^2*/ 
	      zmatind[s] =ind_y_2 + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	      /*y_ii >= 2l_ix_i - l_i^2*/ 
	      zmatind[s] =ind_y_3 + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	    }
	
	  if(j>i)
	    {	  
	      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
	      zmatind[s] = ind_ineg_q + i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=1; 
	      s++;   
	      
	      /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
	      zmatind[s] =ind_y_1+ i* (2*(cqp->n) -i -1)/2 + j-i-1;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
	      zmatind[s] =ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	      /*y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
	      zmatind[s] =ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	    }
	}
 
  

   /* Second member of constrains */
   zrhs = alloc_vector_d (nbrows);
  i=0;
   
  /* Ax = b */
  for(j=0;j<cqp->m;j++)
    {
      zrhs[i]=cqp->b[j];
      i++;
    }
  /* Dx <= e*/
  for(j=0;j<cqp->p;j++)
    {
      zrhs[i]=cqp->e[j];
      i++;
    }
   
  /*Aq,y + Aq_0 x = bq*/
  for(j=0;j<cqp->mq;j++)
    {
      zrhs[i]=cqp->bq[j];
      i++;
    }
  
  /*Dq,y + Dq_0 x <= eq*/  
  for(j=0;j<cqp->pq;j++)
    {
      zrhs[i]=cqp->eq[j];
      i++;
    }

  
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=-(double)cqp->u[k] *(double)cqp->l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->n;j++)
    for(k=j+1;k<cqp->n;k++)
   	{	  
	  zrhs[i]=-(double)cqp->u[j] *(double)cqp->l[k] ;
	  i++;
	}
   
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=(double)cqp->u[k] *(double)cqp->u[j] ;
	  i++;
	}
      
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=(double)cqp->l[k] *(double)cqp->l[j] ;
	  i++;
	}
  
  /* Sense of constraints */
  zsense    = alloc_string(nbrows);
   
  for(i=0;i<ind_eg;i++)
    zsense[i]= 'E';
  for(i;i<ind_ineg;i++)
    zsense[i]= 'L';
  for(i;i<ind_eg_q;i++)
    zsense[i]= 'E';
  for(;i<nbrows;i++)
    zsense[i]= 'L';
 
   
  /***********************************************************************************/
  /*********************************** Bounds*****************************************/
  /***********************************************************************************/
  zlb= alloc_vector_d(nbcols);
  zub = alloc_vector_d(nbcols);
  
  for(i=0;i<nbcols-1;i++)
    {
            
      zlb[i]=(double)cqp->l[i];
      
      zub[i]=(double)cqp->u[i];
    }

  zlb[i]=cqp->cons;
  zub[i]=cqp->cons;

  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
   
   
  /* Be sure that the new matrix is SDP, compute new lambda_min and add it to diagonal terms*/
  
 
  
  zobj = alloc_vector_d (nbcols);
  for (i=0;i<nbcols;i++)
    zobj[i]=0;

  if (sense_opt == 0)
    zobj[ind_var]=1;
  else
    zobj[ind_var]=-1;
      
 
 
  /**************************************************************************************/
  /******************************** Type of variables ***********************************/
  /**************************************************************************************/
   
  zctype = alloc_string(nbcols); 
   
  
  for(i=0;i<nbcols;i++)
    zctype[i]= 'C';
   
  strcpy (zprobname, "obbt_cplex_var");
   
  *numcols_p   = nbcols;
*numrows_p   = nbrows;
  *objsen_p    = CPX_MIN;   /* The problem is minimization */
   
  *probname_p  = zprobname;
  *obj_p       = zobj;
  *rhs_p       = zrhs;
  *sense_p     = zsense;
  *matbeg_p    = zmatbeg;
  *matcnt_p    = zmatcnt;
  *matind_p    = zmatind;
  *matval_p    = zmatval;
  *lb_p        = zlb;
  *ub_p        = zub;
  *ctype_p = zctype; 


  return (status);
   
}  /* END setproblemdata */



int compute_cplex_obbt_var(int ind_var, int sense_opt, double upperbound)
{
/* Declare pointers for the variables and arrays that will contain
      the data which define the LP problem.  The setproblemdata() routine
      allocates space for the problem data.  */
  int is_improved=0;
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
   int      *linnzcnt= NULL;
   int      *quadnzcnt= NULL;
   double   *qrhs= NULL;
   char     *qsense= NULL;
   int      *linind= NULL;
   double   *linval= NULL;
   int      *qconrow = NULL;
   int      *qconcol = NULL;
   double   *qconval = NULL;

   
   double obbt_sol_adm;
   
   int           cur_numrows, cur_numcols,i;
   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat,k;
   double   objval;
   int nbcols = cqp->n + (cqp->n+1)*cqp->n/2+1;
   double * local_sol = alloc_vector_d (nbcols);
   
   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL;
   int           status;

  

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);
   
   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

   if ( env == NULL ) {
   char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return 2;
   }

  
   /* Fill in the data for the problem.  */

   status = setcqpproblemdata_cplex_obbt_var (&probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval, &lb, &ub,&ctype, ind_var,sense_opt);
  
    if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return 2;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, probname);

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of
      failure, an error message will have been written to the error
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create problem.\n");
      return 2;
   }

   /* Now copy the LP part of the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      return 2;
   }

  
   /* /\* Now add quadratic constraints *\/ */


   if (OBBT_QUAD ==1)
     {
   status = setquadconst_cplex_obbt ( upperbound,&linnzcnt, &quadnzcnt,&qrhs, &qsense, &linind, &linval,&qconrow, &qconcol, &qconval);
     if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return 2;
   }
   
   int* cur_linind;
   double * cur_linval;
   int* cur_qconrow;
   int* cur_qconcol;
   double* cur_qconval;
   int cpt_lin=0;
   int cpt_quad=0;
   
   
   cur_linind = &linind[0];
   cur_linval = &linval[0];
         
   cur_qconrow = &qconrow[0] ;
   cur_qconcol = &qconcol[0];
   cur_qconval = &qconval[0];
   
   status = CPXaddqconstr (env, lp,linnzcnt[0],quadnzcnt[0],qrhs[0],qsense[0], cur_linind,cur_linval, cur_qconrow, cur_qconcol, cur_qconval, NULL);
   if ( status )
     {
       fprintf (stderr, "Failed to copy quadratic constraint.\n");
       return 2;
     }
 
   
     }
   
   
   /* Optimize the problem and obtain solution. */

   status = CPXqpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize QP.\n");
      return 2;
   }
   
   
   solstat = CPXgetstat (env, lp);
   
   if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL &&  solstat!=CPX_STAT_FIRSTORDER)
     {
       printf ( "\n\nCompute obbt sol failed, exit ...\n\n solstat : %d \n\n",solstat);
       return 2;
     }
   
  
  status = CPXgetobjval (env, lp, &obbt_sol_adm);
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return 2;
    }

  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  /*status = CPXgetx (env, lp, cqp->x_local, 0, cur_numcols-1);*/
  status = CPXgetx (env, lp, local_sol, 0, cur_numcols-1);
  
  /*for (i=0;i<qp->n;i++)
    x_local[0][i] = qp->local_sol[i];*/
   
  if ( status ) 
    {
      printf ( "Failed to get optimal x.\n");
      return 2;
    }
   
 
 
  
  if (sense_opt==0)
    {
      if ( obbt_sol_adm -  cqp->l[ind_var] >  EPS_OBBT)
	{
	  is_improved=1;
	  cqp->l[ind_var]=obbt_sol_adm;
	  if ((ind_var < cqp->nb_int) &&  (!is_integer(cqp->l[ind_var])))
	    cqp->l[ind_var] = (int)cqp->l[ind_var]+1;
	}
    }
  else
    {
      obbt_sol_adm = -obbt_sol_adm;
      if (  cqp->u[ind_var] - obbt_sol_adm  >  EPS_OBBT)
	{
	  is_improved=1;
	  cqp->u[ind_var]=obbt_sol_adm;
	  if ((ind_var < cqp->nb_int) &&  (!is_integer(cqp->u[ind_var])))
	    cqp->u[ind_var] = (int)cqp->u[ind_var];
	}
    }
  
  

  
 
  /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */

      if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
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
   free_and_null ((char **) &lb);
   free_and_null ((char **) &ub);
   free_and_null ((char **) &ctype);
  
   return is_improved;

}






int setqpproblemdata_cplex_obbt_cont (char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, double **lb_p, double **ub_p, char ** ctype_p, int ind_cont,int sense_opt)
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
   
  int i,j,k,l,s;

  double * new_q;

  zprobname = alloc_string(16); 
  
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
   // + 1 variable for the constant term
  int nbcols = cqp->n + (cqp->n)*(cqp->n+1)/2 +1;
  zmatbeg   = alloc_vector(nbcols);   
  zmatcnt   = alloc_vector(nbcols); 

 
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i<=cqp->n;i++)
    zmatbeg[i]= zmatbeg[i-1] + 4*(cqp->n-1) + 3 + cqp->m+cqp->p +cqp->mq +cqp->pq-1;
    
    
  /* y  */
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
      {
	zmatbeg[i]=zmatbeg[i-1];
	if (j==k)
	  zmatbeg[i]= zmatbeg[i] +3;
	else
	  zmatbeg[i]= zmatbeg[i] + 4;
	
	if (!Zero_BB(cqp->mq))
	  zmatbeg[i]= zmatbeg[i] + cqp->mq;
	if (!Zero_BB(cqp->pq)  && (ind_cont >= cqp->p)) 
	  zmatbeg[i]= zmatbeg[i] + cqp->pq-1;
	else
	  if(!Zero_BB(cqp->pq))
	    zmatbeg[i]= zmatbeg[i] + cqp->pq;
	i++;
      }
      
 
 
  int numnz;
 
  numnz = zmatbeg[nbcols-1];
  
  for(i=0;i<nbcols-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  
  zmatcnt[nbcols - 1] =0;

  
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
   

  /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  int ind_ineg,ind_ineg_q;

 
  int ind_eg = cqp->m;
  if (ind_cont < cqp->p)
    ind_ineg = ind_eg+cqp->p-1;
  else
    ind_ineg = ind_eg+cqp->p;
  
  int ind_eg_q = ind_ineg + cqp->mq;
 
  if (ind_cont >= cqp->p)
    ind_ineg_q=ind_eg_q + cqp->pq-1;
  else
    ind_ineg_q = ind_eg_q + cqp->pq;
  
  int ind_y_1  = ind_ineg_q + (cqp->n+1)*cqp->n/2;
  int ind_y_2 = ind_y_1 + (cqp->n-1)*cqp->n/2 ;
  int ind_y_3 = ind_y_2 + (cqp->n+1)*cqp->n/2 ;
  int nbrows = ind_y_3 + (cqp->n+1)*cqp->n/2 ;

  int ind_cont_k=0;
  /* Contraints on x in order x_1 x_2 ... x_n */
  for(i=0;i<cqp->n;i++)
    {
      
      /*Ax = b*/
      for(j= 0;j<ind_eg;j++)
	{ 
	  
	  zmatval[s]=cqp->a[ij2k(j,i,cqp->n)];
	  zmatind[s]=j;
	  s++;
	}
      ind_cont_k=ind_eg;
      /* Dx =e*/
      for(j=ind_eg;j<cqp->m+cqp->p;j++)
	if (ind_cont!=j)
	    {
	      zmatval[s]=cqp->d[ij2k(j-ind_eg,i,cqp->n)];
	      zmatind[s]=ind_cont_k;
	      ind_cont_k++;
	      s++;
	    }
		
      
      /*Aq,y + Aq_0 x = bq*/
      for(j=ind_ineg;j<ind_eg_q;j++)
	{
	  zmatval[s]=2*cqp->aq[ijk2l(j-ind_ineg,0,i+1,cqp->n+1,cqp->n+1)];
	  zmatind[s]=j;
	  s++;
	}
    
      ind_cont_k=ind_eg_q;
      /*Dq,y + Dq_0 x <= eq*/
      for(j=ind_eg_q;j<cqp->m+cqp->p+cqp->mq+cqp->pq;j++)
	if (j!=ind_cont)
	  {
	    zmatval[s]=2*cqp->dq[ijk2l(j- cqp->m-cqp->p-cqp->mq,0,i+1,cqp->n+1,cqp->n+1)];
	    zmatind[s]=ind_cont_k;
	    ind_cont_k++;
	    s++;
	  }
	
        
      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
      for (j=0;j<cqp->n;j++)
	{
	  if(j>i)
	      { 
		zmatind[s]= ind_ineg_q + i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=-(double)cqp->u[j];
		s++;
	      }
	  if(i>j)
	    { 
		zmatind[s]=ind_ineg_q + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=-(double)cqp->l[j];
		s++;
	    } 
	  if(i==j)
	    {
		zmatind[s]=ind_ineg_q + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=-((double)cqp->l[i]+(double)cqp->u[j]);
		s++;
	      }
	}
      
	        
      /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
      for (j=0;j<cqp->n;j++)
	{
	  if(j>i)
	    { 
	      zmatind[s]= ind_y_1 + i* (2*(cqp->n) - i -1)/2 + j-i-1;
	      zmatval[s]=-(double)cqp->l[j];
	      s++;
	    }
	  if(i>j)
	    { 
	      zmatind[s]=ind_y_1 + j* (2*(cqp->n) -j+ -1)/2  + i-j-1;
	      zmatval[s]=-(double)cqp->u[j];
	      s++;
	    } 
	  }
      
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<cqp->n;j++)
	{
	    if(j>i)
	      { 
		zmatind[s] =ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)cqp->u[j]; 
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=ind_y_2 + j* (2*(cqp->n) -j +1)/2  + i-j;
		zmatval[s]=(double)cqp->u[j];
	      s++;
	      }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_2  + j* (2*(cqp->n) -j +1)/2;
		zmatval[s]=2*(double)cqp->u[i];
		s++;
	      }
	  }
      
      /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
      for (j=0;j<cqp->n;j++)
	  {
	    if(j>i)
	      { 
		zmatind[s] =ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)cqp->l[j]; 
		s++;
	      }
	    if(i>j)
	    { 
	      zmatind[s]=ind_y_3 + j* (2*(cqp->n) -j +1)/2  + i-j;
	      zmatval[s]=(double)cqp->l[j];
	      s++;
	    }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_3  + j* (2*(cqp->n) -j +1)/2;
		zmatval[s]=2*(double)cqp->l[i];
		s++;
	      }
	  }
    }

      
  /* Constraints on y in order y_11 y_12 .. y_1n , ... ,y_n1 .. y_nn */   
  int ind_contk=cqp->m+cqp->p+cqp->mq;
  for(i=0;i<cqp->n;i++)
    for(j=i;j<cqp->n;j++)
      {
	  
	/*Aq,y + Aq_0 x = bq*/
	for(k=ind_ineg;k<ind_eg_q;k++)
	  {
	      if (i==j)
	      {
		zmatval[s]=cqp->aq[ijk2l(k-ind_ineg,i+1,j+1,cqp->n+1,cqp->n+1)];
		zmatind[s]=k;
		s++;
	  	}
	    else
	      {
	  	  zmatval[s]=2*cqp->aq[ijk2l(k-ind_ineg,i+1,j+1,cqp->n+1,cqp->n+1)];
	  	  zmatind[s]=k;
	  	  s++;
	      }
	  }
	ind_cont_k=ind_eg_q;
	/* Dq,y + Dq_0 x <=eq*/
	for(k=ind_eg_q;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
	  if (k != ind_cont)
	    {
	      if (i==j)
		{
		  zmatval[s]=cqp->dq[ijk2l(k-cqp->m-cqp->p-cqp->mq,i+1,j+1,cqp->n+1,cqp->n+1)];
		  zmatind[s]=ind_cont_k;
		  ind_cont_k++;
		  s++;
		  }
	      else
		{
		  zmatval[s]=2*cqp->dq[ijk2l(k-cqp->m-cqp->p-cqp->mq,i+1,j+1,cqp->n+1,cqp->n+1)];
		  zmatind[s]=ind_cont_k;
		  ind_cont_k++;
		  s++;
		}
	      }
	  
 
	  if(i==j)
	    {	  
	      /*y_ii <= l_ix_i+u_ix_i - l_iu_i*/ 
	      zmatind[s] =ind_ineg_q + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ii >= 2u_ix_i - u_i^2*/ 
	      zmatind[s] =ind_y_2 + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	      /*y_ii >= 2l_ix_i - l_i^2*/ 
	      zmatind[s] =ind_y_3 + i* (2*(cqp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	    }
	
	  if(j>i)
	    {	  
	      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
	      zmatind[s] = ind_ineg_q + i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=1; 
	      s++;   
	      
	      /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
	      zmatind[s] =ind_y_1+ i* (2*(cqp->n) -i -1)/2 + j-i-1;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
	      zmatind[s] =ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	      /*y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
	      zmatind[s] =ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	    }
	}
 
  
 
   /* Second member of constrains */
   zrhs = alloc_vector_d (nbrows);
  i=0;
   
  /* Ax = b */
  for(j=0;j<cqp->m;j++)
    {
	zrhs[i]=cqp->b[j];
	i++;
      }
  
  /* Dx <= e*/
  for(j=0;j<cqp->p;j++)
    if(j!=ind_cont-cqp->m)
      {
	zrhs[i]=cqp->e[j];
	i++;
      }
   
  /*Aq,y + Aq_0 x = bq*/
  for(j=0;j<cqp->mq;j++)
    {
	zrhs[i]=cqp->bq[j];
	i++;
      }
  
  /*Dq,y + Dq_0 x <= eq*/  
  for(j=0;j<cqp->pq;j++)
    if (j!=ind_cont-cqp->m-cqp->p-cqp->mq)
      {
	zrhs[i]=cqp->eq[j];
	i++;
      }

  
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=-(double)cqp->u[k] *(double)cqp->l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->n;j++)
    for(k=j+1;k<cqp->n;k++)
   	{	  
	  zrhs[i]=-(double)cqp->u[j] *(double)cqp->l[k] ;
	  i++;
	}
   
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=(double)cqp->u[k] *(double)cqp->u[j] ;
	  i++;
	}
      
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n;k++)
   	{	  
	  zrhs[i]=(double)cqp->l[k] *(double)cqp->l[j] ;
	  i++;
	}
  
  /* Sense of constraints */
  zsense    = alloc_string(nbrows);
   
  for(i=0;i<ind_eg;i++)
    zsense[i]= 'E';
  for(i;i<ind_ineg;i++)
    zsense[i]= 'L';
  for(i;i<ind_eg_q;i++)
    zsense[i]= 'E';
  for(;i<nbrows;i++)
    zsense[i]= 'L';
 
   
  /***********************************************************************************/
  /*********************************** Bounds*****************************************/
  /***********************************************************************************/
  zlb= alloc_vector_d(nbcols);
  zub = alloc_vector_d(nbcols);
  
  for(i=0;i<nbcols-1;i++)
    {
      zlb[i]=(double)cqp->l[i];
      zub[i]=(double)cqp->u[i];
    }
  
  zlb[i]=cqp->cons;
  zub[i]=cqp->cons;
  
  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
   
   
  /* optimize the constraint ind_cont*/

   
  zobj = alloc_vector_d (nbcols);
 
  if (ind_cont < cqp->p)
    {
      for (i=0;i<cqp->n;i++)
	zobj[i]=-cqp->d[ij2k(ind_cont-cqp->m,i,cqp->n)];
      for (i;i<nbcols;i++)
	zobj[i]=0;
    }

  if (ind_cont >= cqp->p)
    {
      for (i=0;i<cqp->n;i++)
	zobj[i]=-2*cqp->dq[ijk2l(ind_cont-cqp->m-cqp->p-cqp->mq,0,i+1,cqp->n+1,cqp->n+1)];
      
      for (j=0;j<cqp->n;j++)
	for (k=j;k<cqp->n;k++)
	  {
	    if (j==k)
	      zobj[i]=-cqp->dq[ijk2l(ind_cont-cqp->m-cqp->p-cqp->mq,j+1,k+1,cqp->n+1,cqp->n+1)];
	    else
	      zobj[i]=-2*cqp->dq[ijk2l(ind_cont-cqp->m-cqp->p-cqp->mq,j+1,k+1,cqp->n+1,cqp->n+1)];
	    i++;
	  }
    }


  
  
 
 
  /**************************************************************************************/
  /******************************** Type of variables ***********************************/
  /**************************************************************************************/
   
  zctype = alloc_string(nbcols); 
   
  
  for(i=0;i<nbcols;i++)
    zctype[i]= 'C';
   
  strcpy (zprobname, "obbt_cplex_cont");
   
  *numcols_p   = nbcols;
*numrows_p   = nbrows;
  *objsen_p    = CPX_MIN;   /* The problem is minimization */
   
  *probname_p  = zprobname;
  *obj_p       = zobj;
  *rhs_p       = zrhs;
  *sense_p     = zsense;
  *matbeg_p    = zmatbeg;
  *matcnt_p    = zmatcnt;
  *matind_p    = zmatind;
  *matval_p    = zmatval;
  *lb_p        = zlb;
  *ub_p        = zub;
  *ctype_p = zctype; 


  return (status);
   
}  /* END setproblemdata */



int compute_cplex_obbt_cont(int ind_cont, int sense_opt,double upperbound)
{
/* Declare pointers for the variables and arrays that will contain
      the data which define the LP problem.  The setproblemdata() routine
      allocates space for the problem data.  */
  int is_improved=0;
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
   int      *linnzcnt= NULL;
   int      *quadnzcnt= NULL;
   double   *qrhs= NULL;
   char     *qsense= NULL;
   int      *linind= NULL;
   double   *linval= NULL;
   int      *qconrow = NULL;
   int      *qconcol = NULL;
   double   *qconval = NULL;
   
   double obbt_sol_adm;
   
   int           cur_numrows, cur_numcols,i;
   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat,k;
   double   objval;
   int nbcols = qp->n + (qp->n+1)*qp->n/2 +1 ;
   double * local_sol = alloc_vector_d (nbcols);
   
   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL;
   int           status;

  

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);
   
   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

   if ( env == NULL ) {
   char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return 2;
   }

   /* Turn on output to the screen */
  

   /* Fill in the data for the problem.  */

   status = setqpproblemdata_cplex_obbt_cont (&probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval, &lb, &ub,&ctype, ind_cont,sense_opt);
  
    if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return 2;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, probname);

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of
      failure, an error message will have been written to the error
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create problem.\n");
      return 2;
   }

   /* Now copy the LP part of the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      return 2;
   }

  

   
   /* /\* Now add quadratic constraints *\/ */
   if (OBBT_QUAD == 1)
     {
   status = setquadconst_cplex_obbt ( upperbound,&linnzcnt, &quadnzcnt,&qrhs, &qsense, &linind, &linval,&qconrow, &qconcol, &qconval);
     if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return 2;
   }
   
   int* cur_linind;
   double * cur_linval;
   int* cur_qconrow;
   int* cur_qconcol;
   double* cur_qconval;
   int cpt_lin=0;
   int cpt_quad=0;
   
   
   cur_linind = &linind[0];
   cur_linval = &linval[0];
         
   cur_qconrow = &qconrow[0] ;
   cur_qconcol = &qconcol[0];
   cur_qconval = &qconval[0];
   
   status = CPXaddqconstr (env, lp,linnzcnt[0],quadnzcnt[0],qrhs[0],qsense[0], cur_linind,cur_linval, cur_qconrow, cur_qconcol, cur_qconval, NULL);
   if ( status )
     {
       fprintf (stderr, "Failed to copy quadratic constraint.\n");
       return 2;
     }
     }
 
   
    
   /* Optimize the problem and obtain solution. */

   status = CPXqpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize QP.\n");
      return 2;
   }
   
  
   solstat = CPXgetstat (env, lp);
   
   if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL &&  solstat!=CPX_STAT_FIRSTORDER)
     {
       printf ( "\n\nCompute obbt with cont %d failed, exit ...\n\n solstat : %d \n\n",ind_cont,solstat);
       return 2;
     }
   
  
  status = CPXgetobjval (env, lp, &obbt_sol_adm);
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return 2;
    }

  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  /*status = CPXgetx (env, lp, cqp->x_local, 0, cur_numcols-1);*/
  status = CPXgetx (env, lp, local_sol, 0, cur_numcols-1);
  
  /*for (i=0;i<qp->n;i++)
    x_local[0][i] = qp->local_sol[i];*/
   
  if ( status ) 
    {
      printf ( "Failed to get optimal x.\n");
      return 2;
    }
   
 
  obbt_sol_adm=-obbt_sol_adm;
  if (ind_cont < cqp->p + cqp->m)
    if ( obbt_sol_adm - cqp->e[ind_cont-cqp->m]  < - EPS_OBBT )
      {
	is_improved=1;
	//cqp->e[ind_cont-cqp->m]=obbt_sol_adm;
      }
  
  if (ind_cont >=cqp->p + cqp->m + cqp->mq)
    if (obbt_sol_adm -  cqp->eq[ind_cont-cqp->m-cqp->p-cqp->mq] < - EPS_OBBT)
      {
	is_improved=1;
	//cqp->eq[ind_cont-cqp->m-cqp->p-cqp->mq]=obbt_sol_adm;
      }
       
   

  /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      
      if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
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
   free_and_null ((char **) &lb);
   free_and_null ((char **) &ub);
   free_and_null ((char **) &ctype);
  
   return is_improved;

}



void presolve_with_obbt(int * nb_pass_var, int * nb_pass_cont,double upperbound, int nb_pass){

  int proceed_obbt=1;
  int improved_inf, improved_sup;
  int nbpass_var=0;
  int improved_cont;
  int nbpass_cont=0;
  int nbcols = cqp->n + (cqp->n)*(cqp->n+1)/2;
  int nb_cont_ineg_lin=cqp->p;
  int nb_cont_ineg_quad=cqp->pq;
  int i;

  

  /*For variables bounds*/

  while(proceed_obbt == 1 && nbpass_var <nb_pass)
    {
      
     
      
      proceed_obbt =0;
      improved_inf=0;
      improved_sup=0;
      for(i=0;i<nbcols;i++)
	{
	  improved_inf=compute_cplex_obbt_var(i,0,upperbound);
	  if (improved_inf ==1)
	    {
	      proceed_obbt =1;
	    
	    }
	}
      for(i=0;i<cqp->n;i++)
	{
	  improved_sup=compute_cplex_obbt_var(i,1,upperbound);
	  if (improved_sup ==1)
	    {
	      proceed_obbt =1;
	      
	    }
	}
      nbpass_var++;

     


      /* for(i=0;i<cqp->nb_int;i++) */
      /* 	{ */
      /* 	  if (!is_integer(cqp->u[i])) */
      /* 	    cqp->u[i] = (int)cqp->u[i]; */
      /* 	  if (!is_integer(cqp->l[i])) */
      /* 	    cqp->l[i] = (int)cqp->l[i] +1; */
      /* 	} */
     
      
      improved_cont=0;
      if (OBBT_CONT ==1)
	{
	  for(i=cqp->m;i<nb_cont_ineg_lin+cqp->m;i++)
	    {
	      improved_cont=compute_cplex_obbt_cont(i,1,upperbound);
	      if (improved_cont ==1)
		{
		  proceed_obbt =1;
		  reduce_c_miqcp_obbt( i);
		  i--;
		  nb_cont_ineg_lin = cqp->p;
		}
	    }
	  
	  for(i=cqp->m+cqp->p+cqp->mq;i<cqp->m+cqp->p+cqp->mq+nb_cont_ineg_quad;i++)
	    {
	      improved_cont=compute_cplex_obbt_cont(i,1,upperbound);
	      if (improved_cont ==1)
		{
		  proceed_obbt =1;
		  reduce_c_miqcp_obbt( i);
		  i--;
		  nb_cont_ineg_quad = cqp->pq;
		}
	    }
	  nbpass_cont++;
	}
      *nb_pass_var=nbpass_var;
      *nb_pass_cont=nbpass_cont;
    }
}




int setcqpproblemdata_cplex_obbt_init_bounds ( char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, double **lb_p, double **ub_p, char ** ctype_p, int ind_var,int sense_opt)
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
   
  int i,j,k,l,s;

  double * new_q;

  zprobname = alloc_string(16); 
  
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
  // + 1 variable for the constant term
  int nbcols = qp->n + (qp->n)*(qp->n+1)/2 + 1;
  zmatbeg   = alloc_vector(nbcols);   
  zmatcnt   = alloc_vector(nbcols); 

 
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i<=qp->n;i++)
    zmatbeg[i]= zmatbeg[i-1] + 4*(qp->n-1) + 3 + qp->m+qp->p +qp->mq +qp->pq;
    
    
  /* y  */
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
      {
	zmatbeg[i]=zmatbeg[i-1];
	if (j==k)
	  zmatbeg[i]= zmatbeg[i] +3;
	else
	  zmatbeg[i]= zmatbeg[i] + 4;
	
	if (!Zero_BB(qp->mq))
	  zmatbeg[i]= zmatbeg[i] + qp->mq;
	if (!Zero_BB(qp->pq))
	  zmatbeg[i]= zmatbeg[i] + qp->pq;
	i++;
      }
      
 
 
  int numnz;
  
  numnz = zmatbeg[nbcols-1];
  
  for(i=0;i<nbcols-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  zmatcnt[nbcols - 1] = 0;
  

  
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
  
  /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  int ind_eg = qp->m;
  int ind_ineg = ind_eg+qp->p;
  int ind_eg_q = ind_ineg + qp->mq;
  int ind_ineg_q = ind_eg_q + qp->pq;
  int ind_y_1  = ind_ineg_q + (qp->n+1)*qp->n/2;
  int ind_y_2 = ind_y_1 + (qp->n-1)*qp->n/2 ;
  int ind_y_3 = ind_y_2 + (qp->n+1)*qp->n/2 ;
  int nbrows = ind_y_3 + (qp->n+1)*qp->n/2 ;
 
  /* Contraints on x in order x_1 x_2 ... x_n */
  for(i=0;i<qp->n;i++)
    {
      
      /*Ax = b*/
      for(j= 0;j<ind_eg;j++)
	{ 
	  zmatval[s]=qp->a[ij2k(j,i,qp->n)];
	  zmatind[s]=j;
	  s++;
	}
       
      /* Dx =e*/
      for(j=ind_eg;j<ind_ineg;j++)
	{
	  zmatval[s]=qp->d[ij2k(j-ind_eg,i,qp->n)];
	  zmatind[s]=j;
	  s++;
	}
      
      /*Aq,y + Aq_0 x = bq*/
      for(j=ind_ineg;j<ind_eg_q;j++)
      	{
      	  zmatval[s]=2*qp->aq[ijk2l(j-ind_ineg,0,i+1,qp->n+1,qp->n+1)];
      	  zmatind[s]=j;
      	  s++;
      	}
      
      /*Dq,y + Dq_0 x <= eq*/
      for(j=ind_eg_q;j<ind_ineg_q;j++)
	{
	  zmatval[s]=2*qp->dq[ijk2l(j- ind_eg_q,0,i+1,qp->n+1,qp->n+1)];
	  zmatind[s]=j;
	  s++;
	} 
        
      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
      for (j=0;j<qp->n;j++)
	{
	  if(j>i)
	      { 
		zmatind[s]= ind_ineg_q + i* (2*(qp->n) -i +1)/2 + j-i;
		zmatval[s]=-(double)qp->u[j];
		s++;
	      }
	  if(i>j)
	    { 
		zmatind[s]=ind_ineg_q + j* (2*(qp->n) -j +1)/2  + i-j;
		zmatval[s]=-(double)qp->l[j];
		s++;
	    } 
	  if(i==j)
	    {
		zmatind[s]=ind_ineg_q + j* (2*(qp->n) -j +1)/2  + i-j;
		zmatval[s]=-((double)qp->l[i]+(double)qp->u[j]);
		s++;
	      }
	}
      
	        
      /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
      for (j=0;j<qp->n;j++)
	{
	  if(j>i)
	    { 
	      zmatind[s]= ind_y_1 + i* (2*(qp->n) - i -1)/2 + j-i-1;
	      zmatval[s]=-(double)qp->l[j];
	      s++;
	    }
	  if(i>j)
	    { 
	      zmatind[s]=ind_y_1 + j* (2*(qp->n) -j+ -1)/2  + i-j-1;
	      zmatval[s]=-(double)qp->u[j];
	      s++;
	    } 
	  }
      
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<qp->n;j++)
	{
	    if(j>i)
	      { 
		zmatind[s] =ind_y_2+ i* (2*(qp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)qp->u[j]; 
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=ind_y_2 + j* (2*(qp->n) -j +1)/2  + i-j;
		zmatval[s]=(double)qp->u[j];
	      s++;
	      }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_2  + j* (2*(qp->n) -j +1)/2;
		zmatval[s]=2*(double)qp->u[i];
		s++;
	      }
	  }
      
      /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
      for (j=0;j<qp->n;j++)
	  {
	    if(j>i)
	      { 
		zmatind[s] =ind_y_3+ i* (2*(qp->n) -i +1)/2 + j-i;
		zmatval[s]=(double)qp->l[j]; 
		s++;
	      }
	    if(i>j)
	    { 
	      zmatind[s]=ind_y_3 + j* (2*(qp->n) -j +1)/2  + i-j;
	      zmatval[s]=(double)qp->l[j];
	      s++;
	    }   
	    if(i==j)
	      {
		zmatind[s]=ind_y_3  + j* (2*(qp->n) -j +1)/2;
		zmatval[s]=2*(double)qp->l[i];
		s++;
	      }
	  }
    }

      
  /* Constraints on y in order y_11 y_12 .. y_1n , ... ,y_n1 .. y_nn */   
  
  for(i=0;i<qp->n;i++)
    for(j=i;j<qp->n;j++)
      {
	  
	/*Aq,y + Aq_0 x = bq*/
	for(k=ind_ineg;k<ind_eg_q;k++)
	  {
	    if (i==j)
	      {
		zmatval[s]=qp->aq[ijk2l(k-ind_ineg,i+1,j+1,qp->n+1,qp->n+1)];
		zmatind[s]=k;
		s++;
	  	}
	    else
	      {
	  	  zmatval[s]=2*qp->aq[ijk2l(k-ind_ineg,i+1,j+1,qp->n+1,qp->n+1)];
	  	  zmatind[s]=k;
	  	  s++;
	      }
	  }
	  
	  /* Dq,y + Dq_0 x <=eq*/
	  for(k=ind_eg_q;k<ind_ineg_q;k++)
	    {
	      if (i==j)
		{
		  zmatval[s]=qp->dq[ijk2l(k-ind_eg_q,i+1,j+1,qp->n+1,qp->n+1)];
		  zmatind[s]=k;
		  s++;
		}
	      else
		{
		  zmatval[s]=2*qp->dq[ijk2l(k-ind_eg_q,i+1,j+1,qp->n+1,qp->n+1)];
		  zmatind[s]=k;
		  s++;
		}
	    }

 
	  if(i==j)
	    {	  
	      /*y_ii <= l_ix_i+u_ix_i - l_iu_i*/ 
	      zmatind[s] =ind_ineg_q + i* (2*(qp->n) -i +1)/2;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ii >= 2u_ix_i - u_i^2*/ 
	      zmatind[s] =ind_y_2 + i* (2*(qp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	      /*y_ii >= 2l_ix_i - l_i^2*/ 
	      zmatind[s] =ind_y_3 + i* (2*(qp->n) -i +1)/2;
	      zmatval[s]=-1; 
	      s++;   
	      
	    }
	
	  if(j>i)
	    {	  
	      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
	      zmatind[s] = ind_ineg_q + i* (2*(qp->n) -i +1)/2 + j-i;
	      zmatval[s]=1; 
	      s++;   
	      
	      /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
	      zmatind[s] =ind_y_1+ i* (2*(qp->n) -i -1)/2 + j-i-1;
	      zmatval[s]=1; 
	      s++;   
	      
	      /*y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
	      zmatind[s] =ind_y_2+ i* (2*(qp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	      /*y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
	      zmatind[s] =ind_y_3+ i* (2*(qp->n) -i +1)/2 + j-i;
	      zmatval[s]=-1; 
	      s++;
	      
	    }
	}
 
  

   /* Second member of constrains */
   zrhs = alloc_vector_d (nbrows);
  i=0;
   
  /* Ax = b */
  for(j=0;j<qp->m;j++)
    {
      zrhs[i]=qp->b[j];
      i++;
    }
  /* Dx <= e*/
  for(j=0;j<qp->p;j++)
    {
      zrhs[i]=qp->e[j];
      i++;
    }
   
  /*Aq,y + Aq_0 x = bq*/
  for(j=0;j<qp->mq;j++)
    {
      zrhs[i]=qp->bq[j];
      i++;
    }
  
  /*Dq,y + Dq_0 x <= eq*/  
  for(j=0;j<qp->pq;j++)
    {
      zrhs[i]=qp->eq[j];
      i++;
    }

  
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
   	{	  
	  zrhs[i]=-(double)qp->u[k] *(double)qp->l[j] ;
	  i++;
	}
  
  for(j=0;j<qp->n;j++)
    for(k=j+1;k<qp->n;k++)
   	{	  
	  zrhs[i]=-(double)qp->u[j] *(double)qp->l[k] ;
	  i++;
	}
   
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
   	{	  
	  zrhs[i]=(double)qp->u[k] *(double)qp->u[j] ;
	  i++;
	}
      
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
   	{	  
	  zrhs[i]=(double)qp->l[k] *(double)qp->l[j] ;
	  i++;
	}
  
  /* Sense of constraints */
  zsense    = alloc_string(nbrows);
   
  for(i=0;i<ind_eg;i++)
    zsense[i]= 'E';
  for(i;i<ind_ineg;i++)
    zsense[i]= 'L';
  for(i;i<ind_eg_q;i++)
    zsense[i]= 'E';
  for(;i<nbrows;i++)
    zsense[i]= 'L';
 
   
  /***********************************************************************************/
  /*********************************** Bounds*****************************************/
  /***********************************************************************************/
  zlb= alloc_vector_d(nbcols);
  zub = alloc_vector_d(nbcols);
  
  for(i=0;i<nbcols-1;i++)
    {
            
      zlb[i]=(double)qp->l[i];
      
      zub[i]=(double)qp->u[i];
    }

  zlb[i]=qp->cons;
  zub[i]=qp->cons;

  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
   
   
  /* Be sure that the new matrix is SDP, compute new lambda_min and add it to diagonal terms*/
  
 
  
  zobj = alloc_vector_d (nbcols);
  for (i=0;i<nbcols;i++)
    zobj[i]=0;

  if (sense_opt == 0)
    zobj[ind_var]=1;
  else
    zobj[ind_var]=-1;
      
 
 
  /**************************************************************************************/
  /******************************** Type of variables ***********************************/
  /**************************************************************************************/
   
  zctype = alloc_string(nbcols); 
   
  
  for(i=0;i<nbcols;i++)
    zctype[i]= 'C';
   
  strcpy (zprobname, "obbt_cplex_var");
   
  *numcols_p   = nbcols;
*numrows_p   = nbrows;
  *objsen_p    = CPX_MIN;   /* The problem is minimization */
   
  *probname_p  = zprobname;
  *obj_p       = zobj;
  *rhs_p       = zrhs;
  *sense_p     = zsense;
  *matbeg_p    = zmatbeg;
  *matcnt_p    = zmatcnt;
  *matind_p    = zmatind;
  *matval_p    = zmatval;
  *lb_p        = zlb;
  *ub_p        = zub;
  *ctype_p = zctype; 


  return (status);
   
}  
void compute_cplex_obbt_init_bounds(int ind_var, int sense_opt)
{
/* Declare pointers for the variables and arrays that will contain
      the data which define the LP problem.  The setproblemdata() routine
      allocates space for the problem data.  */
  int is_improved=0;
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


   
   double obbt_sol_adm;
   
   int           cur_numrows, cur_numcols,i;
   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat,k;
   double   objval;
   int nbcols = qp->n + (qp->n)*(qp->n+1)/2 + 1;
   double * local_sol = alloc_vector_d (nbcols);
   
   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL;
   int           status;

  

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);
   
   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

   if ( env == NULL ) {
   char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return ;
   }

  
   /* Fill in the data for the problem.  */

   status = setcqpproblemdata_cplex_obbt_init_bounds (&probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval, &lb, &ub,&ctype, ind_var,sense_opt);
  
    if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return ;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, probname);

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of
      failure, an error message will have been written to the error
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create problem.\n");
      return ;
   }

   /* Now copy the LP part of the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      return ;
   }

     
   
   /* Optimize the problem and obtain solution. */

   status = CPXqpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize QP.\n");
      return ;
   }
   
   
   solstat = CPXgetstat (env, lp);
   
   if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL && solstat!=CPX_STAT_FIRSTORDER)
     {
       printf ( "\n\nCompute obbt sol failed, exit ...\n\n solstat : %d \n\n",solstat);
       return ;
     }
   
  
  status = CPXgetobjval (env, lp, &obbt_sol_adm);
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return ;
    }

  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  /*status = CPXgetx (env, lp, qp->x_local, 0, cur_numcols-1);*/
  status = CPXgetx (env, lp, local_sol, 0, cur_numcols-1);
  
  /*for (i=0;i<qp->n;i++)
    x_local[0][i] = qp->local_sol[i];*/
   
  if ( status ) 
    {
      printf ( "Failed to get optimal x.\n");
      return ;
    }
   
 
 
  
  if (sense_opt==0)
    {
      if ( obbt_sol_adm -  qp->l[ind_var] >  EPS_OBBT)
	{
	  is_improved=1;
	  qp->l[ind_var]=obbt_sol_adm;
	  if ((ind_var < qp->nb_int) &&  (!is_integer(qp->l[ind_var])))
	    qp->l[ind_var] = (int)qp->l[ind_var]+1;
	}
    }
  else
    {
      obbt_sol_adm = -obbt_sol_adm;
      if (  qp->u[ind_var] - obbt_sol_adm  >  EPS_OBBT)
	{
	  is_improved=1;
	  qp->u[ind_var]=obbt_sol_adm;
	  if ((ind_var < qp->nb_int) &&  (!is_integer(qp->u[ind_var])))
	    qp->u[ind_var] = (int)qp->u[ind_var];
	}
    }
  
  

  
 
  /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */

      if ( status ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
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
   free_and_null ((char **) &lb);
   free_and_null ((char **) &ub);
   free_and_null ((char **) &ctype);
  
   

}






void initialize_bounds_with_obbt(){
  
  int nbcols = qp->n;
  int i;
  /*For variables bounds*/
  int improved_inf,improved_sup;
  for(i=0;i<nbcols;i++)
    compute_cplex_obbt_init_bounds(i,0);
    
  for(i=0;i<qp->n;i++)
    compute_cplex_obbt_init_bounds(i,1);

  
  
}




/****************************************************************************************************/
/************************************ OBBT for BB ***************************************************/
/****************************************************************************************************/


void presolve_with_obbt_bb(struct miqcp_bab bab,int * nb_pass_var, int * nb_pass_cont, double upperbound,int nb_pass){

  int proceed_obbt=1;
  int improved_inf, improved_sup;
  int nbpass_var=0;
 int nbpass_cont=0;

  int i;
  /*For variables bounds*/
  while(proceed_obbt == 1 && nbpass_var <nb_pass)
    {
      proceed_obbt =0;
      improved_inf=0;
      improved_sup=0;
      for(i=0;i<cqp->n;i++)
	{
	  improved_inf=compute_cplex_obbt_var(i,0,upperbound);
	  if (improved_inf ==1)
	    {
	      proceed_obbt =1;
	      bab.l[i] = cqp->l[i];
	    }
	}
      for(i=0;i<cqp->n;i++)
	{
	  improved_sup=compute_cplex_obbt_var(i,1,upperbound);
	  if (improved_sup ==1)
	    {
	      proceed_obbt =1;
	      bab.u[i] = cqp->u[i];
	    }
	}
      nbpass_var++;
    }
  
  for(i=0;i<cqp->nb_int;i++)
    {
      if (!is_integer(bab.u[i]))
	bab.u[i] = (int)bab.u[i];
      if (!is_integer(bab.l[i]))
	bab.l[i] = (int)bab.l[i] +1;
    }
  
  
  
  
  
  *nb_pass_var=nbpass_var;
  *nb_pass_cont=nbpass_cont;
}


void branch_and_reduce_bb(struct miqcp_bab bab,int * nb_pass, double * slack, double *pi,double upper_bound, double sol_admissible, int numrows){
  
  int i,j;
  double tmp;
  for(i=0;i<numrows;i++)
    {
      if (slack[i]==0)
	{
	  tmp=1;
	}

    }

}
