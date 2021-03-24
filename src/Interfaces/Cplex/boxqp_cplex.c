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


#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include<string.h>



#include "utilities.h"
#include"quad_prog.h"
#include "boxqp_cplex.h"


extern C_MIQCP cqp;
extern long nodecount;
extern double * best_x;
extern double ub_bb;

/* parameters for CPLEX */
extern int MEMORY;
extern double ABS_GAP;
extern  double REL_GAP;
extern int VARSEL;
extern int TIME_LIMIT_CPLEX;
extern double EPS_BETA;


int setqpproblemdata_boxqp_cplex(char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p,double **lb_p, double **ub_p, char ** ctype_p)
{
  char     *zprobname = NULL;     /* Problem name <= 16 characters */        
  double   *zobj = NULL;
  double   *zrhs = NULL;
  char     *zsense = NULL;
  int      *zmatbeg = NULL;
  int      *zmatcnt = NULL;
  int      *zmatind = NULL;
  double   *zmatval = NULL;
  int      *zqmatbeg = NULL;
  int      *zqmatcnt = NULL;
  int      *zqmatind = NULL;
  double   *zqmatval = NULL;
  
  double   *zlb = NULL;
  double   *zub = NULL;
  int      status = 0;
  char *zctype = NULL;
   
  int i,j,k,l,s;

  double * new_q;

  zprobname = alloc_string(16); 
  
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
  
  zmatbeg   = alloc_vector(cqp->nb_col);   
  zmatcnt   = alloc_vector(cqp->nb_col); 

 
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i<=cqp->ind_max_x;i++)
    zmatbeg[i]= zmatbeg[i-1] + cqp->x_cont[i-1];
    
  /* y  */
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{
	  zmatbeg[i]= zmatbeg[i-1];
	  if (j==k)
	    {
	      zmatbeg[i]= zmatbeg[i] +1;
	      i++;
	    }
	  else
	    {
	      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
		{
		  zmatbeg[i]= zmatbeg[i] +2;
		  i++;
		}
	      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
		{
		  zmatbeg[i]= zmatbeg[i] +2;
		  i++;
		}
	    }
	}
  
 
  int numnz;
  
  numnz = zmatbeg[cqp->nb_col-1];
  
  for(i=0;i<cqp->nb_col-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  

 
  zmatcnt[cqp->nb_col - 1] = 0;
  

  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
 
  /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
   
  /* Contraints on x in order x_1 x_2 ... x_n */
  for(i=0;i<cqp->ind_max_x;i++)
    {
      
      /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)])) || (Negative_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	    if(j>i)
	      { 
		zmatind[s]= cqp->cont[ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		zmatval[s]=-1;
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[ j* (2*(cqp->n) -j -1)/2  + i-j-1];
		zmatval[s]=0;
		s++;
	      } 
	  }
    
	        
      /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))|| (Negative_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	  if(j>i)
	    { 
	      zmatind[s]= cqp->cont[cqp->ind_y_1 + i* (2*(cqp->n) - i -1)/2 + j-i-1];
	      zmatval[s]=0;
	      s++;
	    }
	  if(i>j)
	    { 
	      zmatind[s]=cqp->cont[cqp->ind_y_1 + j* (2*(cqp->n) -j+ -1)/2  + i-j-1];
	      zmatval[s]=-1;
	      s++;
	    } 
	  }
      
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))|| (Positive_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	    if(j>i)
	      { 
		zmatind[s] =cqp->cont[cqp->ind_y_2+ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		zmatval[s]=1;
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[cqp->ind_y_2 + j* (2*(cqp->n) -j -1)/2  + i-j-1];
		zmatval[s]=1;
	      s++;
	      }   
	  }
      
      /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))|| (Positive_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	    if(j>i)
	      { 
		zmatind[s] =cqp->cont[cqp->ind_y_3+ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		zmatval[s]=0;
		s++;
	      }
	    if(i>j)
	    { 
	      zmatind[s]=cqp->cont[cqp->ind_y_3 + j* (2*(cqp->n) -j -1)/2  + i-j-1];
	      zmatval[s]=0;
	      s++;
	    }   
	  }
      
      /* y_ii == x_i*/
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]))
	{
	  zmatind[s]=cqp->cont[cqp->ind_y_4 + i];
	  zmatval[s]=1;
	  s++;
	}
    }

      
  /* Constraints on y in order y_11 y_12 .. y_1n , ... ,y_n1 .. y_nn */   
  
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	{
	  
	  if(i==j)
	    {	  
	      /*y_ii = x_i*/ 
	      zmatind[s] =cqp->cont[cqp->ind_y_4+i];
	      zmatval[s]=-1; 
	      s++;
	    }
	  
	  if(j>i)
	    {	  
	     if (Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		{
		  /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
		  zmatind[s] = cqp->cont[ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		  zmatval[s]=1; 
		  s++;   
		  
		  /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
		  zmatind[s] =cqp->cont[cqp->ind_y_1+ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		  zmatval[s]=1; 
		  s++;   
		}
	       if (Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		{
		  /*y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
		  zmatind[s] =cqp->cont[cqp->ind_y_2+ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		  zmatval[s]=-1; 
		  s++;
		  
		  /*y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
		  zmatind[s] =cqp->cont[cqp->ind_y_3+ i* (2*(cqp->n) -i -1)/2 + j-i-1];
		  zmatval[s]=-1; 
		  s++;
		}
	    }
	}

 
   /* Second member of constrains */
   zrhs = alloc_vector_d (cqp->nb_row);
  i=0;
  
 
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=0 ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=0;
	  i++;
	}
   
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=1 ;
	  i++;
	}
      
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=0 ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    if (!Zero_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]))
      {	  
	if (i<cqp->nb_row)
	  {
	    zrhs[i]=0;
	    i++;
	  }
      }
  
  
  /* Sense of constraints */
  zsense    = alloc_string(cqp->nb_row);
  for(i=0;i<cqp->nb_row - cqp->nb_y_5;i++)
    zsense[i]= 'L';
  for(;i<cqp->nb_row;i++)
    zsense[i]= 'E';
   
  /***********************************************************************************/
  /*********************************** Bounds*****************************************/
  /***********************************************************************************/
 
  zlb= alloc_vector_d(cqp->nb_col);
  for(i=0;i<cqp->nb_col-1;i++)
    zlb[i]=0;
 
  zlb[cqp->nb_col-1] =  cqp->cons; 
  
  zub = alloc_vector_d(cqp->nb_col);
  for(i=0;i<cqp->nb_col-1;i++)
    zub[i]=1;
  
  zub[cqp->nb_col-1] = cqp->cons; 
    
  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
   
   
  /* Be sure that the new matrix is SDP, compute new lambda_min and add it to diagonal terms*/
  
 
  zobj = alloc_vector_d (cqp->nb_col);
  for(i=0;i<cqp->n;i++)
    {
      zobj[i]= cqp->c[i] ;
    }
    
  i=cqp->n;
  for(j=0;j<cqp->n;j++)
    for(k=j;k<cqp->n ;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{
	    if (j==k) 
	      {
		zobj[i]=- cqp->beta[ij2k(j,k,cqp->n)] + cqp->new_lambda_min; 
		i++;
	      }
	    else
	      {
		zobj[i]= - 2*(cqp->beta[ij2k(j,k,cqp->n)]);
		i++;
	      }
	  }
  
  zobj[cqp->nb_col-1]=1;

  int numqnz = nb_non_zero_matrix(cqp->new_q, cqp->n,cqp->n);

  zqmatbeg  = alloc_vector(cqp->nb_col);
  zqmatbeg[0]=0;
  for (i=1;i<=cqp->ind_max_x;i++)
    zqmatbeg[i]=zqmatbeg[i-1]+ cqp->ind_max_x;
  
  for (i;i<cqp->nb_col;i++)
    zqmatbeg[i]=zqmatbeg[cqp->ind_max_x];
   
  zqmatcnt   = alloc_vector(cqp->nb_col);
  for (i=0;i<cqp->ind_max_x;i++)
    zqmatcnt[i] = zqmatbeg[i+1] - zqmatbeg[i];
   
  for (i;i<cqp->nb_col;i++)
    zqmatcnt[i]=0; 

  zqmatind   = alloc_vector(cqp->ind_max_x*cqp->ind_max_x);
  zqmatval   = alloc_vector_d(cqp->ind_max_x*cqp->ind_max_x);
   
  k=0;
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=0;j<cqp->ind_max_x;j++)
      { 
	if (i==j)
	  zqmatval[k] = cqp->new_q[ij2k(i,j,cqp->n)] - ((cqp->new_lambda_min)*2) ;
	else
	  zqmatval[k] = cqp->new_q[ij2k(i,j,cqp->n)];
	zqmatind[k] = j;
	k++;
      }
  
 
 
  /**************************************************************************************/
  /******************************** Type of variables ***********************************/
  /**************************************************************************************/
   
  zctype = alloc_string(cqp->nb_col); 
   
  for(i=0;i<cqp->nb_int;i++)
    zctype[i]= 'B';
   
  for(i;i<cqp->nb_col;i++)
    zctype[i]= 'C';
   
  strcpy (zprobname, "boxqp_cplex");
   
  *numcols_p   = cqp->nb_col;
  *numrows_p   = cqp->nb_row;
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
   
}  /* END setproblemdata */



void boxqp_cplex()
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
   int      *qmatbeg = NULL;
   int      *qmatcnt = NULL;
   int      *qmatind = NULL;
   double   *qmatval = NULL;

  double   *lb = NULL;
  double   *ub = NULL;
  char     *ctype = NULL;
  int      solstat;
  
   
  

  CPXENVptr    env = NULL;
  CPXLPptr      lp = NULL;
  int           status;
   
  int           cur_numrows, cur_numcols;
  int file_sol;
  int stdoutsave;


   /* Initialize the CPLEX environment */
  int k,l,temp;
     
  int *i =  alloc_vector(1);
  int *j = alloc_vector(1);

  
  int test_time = time(NULL);
  
   
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
      
  /* Turn on output to the screen */

  status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  if ( status ) 
    {
      printf ("Failure to turn on screen indicator, error %d.\n", status);
      return;
    }
   

  status = CPXsetintparam (env, 4012, 0);
  if ( status )
    {
      printf ("Failure to change param 4012, error %d.\n", status);
      return;
    }

  status = CPXsetdblparam (env, CPX_PARAM_EPGAP,REL_GAP);
  if ( status ) 
    {
      printf ("Failure to change relative gap, error %d.\n", status);
      return;
    }

  

  status = CPXsetdblparam (env, CPX_PARAM_TILIM,TIME_LIMIT_CPLEX);
  if ( status ) 
    {
      printf ("Failure to change time limit, error %d.\n", status);
      return;
    }
  
  /* Set memory parameter*/
  status = CPXsetdblparam (env, CPX_PARAM_WORKMEM, MEMORY);
  if ( status )
    {
      printf ("Failure to change memory, error %d.\n", status);
      return;
    }
  
  /*var selection strategie*/
  status = CPXsetintparam (env, CPX_PARAM_VARSEL, VARSEL);
  if ( status ) 
    {
      printf ("Failure to change var selection strategie, error %d.\n", status);
      return;
    }

  /* Fill in the data for the problem.  */

  status = setqpproblemdata_boxqp_cplex (&probname, &numcols, &numrows, &objsen, &obj, 
				     &rhs, &sense, &matbeg, &matcnt, &matind, &matval,
				      &qmatbeg, &qmatcnt, &qmatind, &qmatval,
				      &lb, &ub, &ctype);
   
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
   status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
   if ( status ) {
      fprintf (stderr, "Failed to copy quadratic matrix.\n");
      return;
   }

  /* Now copy the ctype array */
  status = CPXcopyctype (env, lp, ctype);
  if ( status )
    {
      printf ( "Failed to copy ctype\n");
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
   
  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  status = CPXgetmipx (env, lp, best_x, 0, cur_numcols-1);
  if ( status ) 
    {
      printf ( "Failed to get optimal integer x.\n");
      return;
    }
  
  
  nodecount = CPXgetnodecnt (env, lp); 
  ub_bb = ub_bb + cqp->cons;
 
  
  
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
  free_and_null ((char **) &lb);
  free_and_null ((char **) &ub);
  free_and_null ((char **) &ctype);

  
}

  


