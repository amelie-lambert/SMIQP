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
#include<fcntl.h>
#include<unistd.h>
#include<string.h>
#include<time.h>

#include "utilities.h"
#include"quad_prog.h"
#include "local_sol_cplex.h"
#include"liste.h"
#include"sbb_boxqp_cplex.h"
#ifdef HAVE_SCIP
#include "local_sol_scip.h"
#endif
#include "local_sol_ipopt.h"

extern C_MIQCP cqp;
extern Liste liste;

extern long start_time;
extern double ub_bb;
extern  double lb_bb;
extern double * best_x;
extern long nodecount;
extern int father;
extern double EPS_BETA;
extern double EPS_BB;
extern int TIME_LIMIT_BB;
extern double MIN_SOL_BB;
extern double MAX_SOL_BB;
extern double EPS_ABS_GAP;
extern int PRINT_LEVEL_NODE;
/* parameters for CPLEX */
extern int MEMORY;
extern double ABS_GAP;
extern  double REL_GAP;




int setqpproblemdata_sbb_boxqp_cplex(struct miqcp_bab bab, char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p, double **lb_p, double **ub_p, char ** ctype_p)
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
	  if (j==k && j<cqp->nb_int)
	    {
	      zmatbeg[i]= zmatbeg[i] +1;
	      i++;
	    }
	  else
	    {
	      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
		zmatbeg[i]= zmatbeg[i] +2;
	      
	      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
		zmatbeg[i]= zmatbeg[i] + 2;
	      
	      i++;
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
		zmatind[s]= cqp->cont[ i* (2*(cqp->n) -i +1)/2 + j-i];
		zmatval[s]=-(double)bab.u[j];
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[ j* (2*(cqp->n) -j +1)/2  + i-j];
		zmatval[s]=-(double)bab.l[j];
		s++;
	      } 
	    if(i==j && i >= cqp->nb_int) 
	      {
		zmatind[s]=cqp->cont[ j* (2*(cqp->n) -j +1)/2  + i-j];
		zmatval[s]=-((double)bab.u[i] + bab.l[i]);
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
		zmatval[s]=-(double)bab.l[j];
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[cqp->ind_y_1 + j* (2*(cqp->n) -j+ -1)/2  + i-j-1];
		zmatval[s]=-(double)bab.u[j];
		s++;
	      } 
	  }
	
      /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))|| (Positive_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	    if(j>i)
	      { 
		zmatind[s] =cqp->cont[cqp->ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i];
		zmatval[s]=(double)bab.u[j]; 
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[cqp->ind_y_2 + j* (2*(cqp->n) -j +1)/2  + i-j];
		zmatval[s]=(double)bab.u[j];
		s++;
	      }   
	    if(i==j&& i >= cqp->nb_int) 
	      {
		zmatind[s]=cqp->cont[cqp->ind_y_2  + j* (2*(cqp->n) -j +1)/2];
		zmatval[s]=2*(double)bab.u[i];
		s++;
	      }
	  }
	
      /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
      for (j=0;j<cqp->ind_max_x;j++)
	if ((Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))|| (Positive_BB(cqp->nb_var_y[ij2k(j,i,cqp->n)])))
	  {
	    if(j>i)
	      { 
		zmatind[s] =cqp->cont[cqp->ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i];
		zmatval[s]=(double)bab.l[j]; 
		s++;
	      }
	    if(i>j)
	      { 
		zmatind[s]=cqp->cont[cqp->ind_y_3 + j* (2*(cqp->n) -j +1)/2  + i-j];
		zmatval[s]=(double)bab.l[j];
		s++;
	      }   
	    if(i==j&& i >= cqp->nb_int) 
	      {
		zmatind[s]=cqp->cont[cqp->ind_y_3  + j* (2*(cqp->n) -j +1)/2];
		zmatval[s]=2*(double)bab.l[i];
		s++;
	      }
	  }
	
      /* y_ii <= x_i*/
      if ((Negative_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]) && (i >= cqp->nb_int)) || (!Zero_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]) && (i < cqp->nb_int)))
	{
	 
	  zmatind[s]=cqp->cont[cqp->ind_y_4 + i];
	  zmatval[s]=-1;
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
	      if (Negative_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]) && i >= cqp->nb_int)
		{
		  /*y_ii <= l_ix_i+u_ix_i - l_iu_i*/ 
		  zmatind[s] = cqp->cont[ i* (2*(cqp->n) -i +1)/2];
		  zmatval[s]=1; 
		  s++;
		}

	      /*y_ii <= x_i*/ 
	      if ( (Negative_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]) && (i >= cqp->nb_int)) || (!Zero_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)]) && (i < cqp->nb_int))) 
		{
		  zmatind[s] =cqp->cont[cqp->ind_y_4+i];
		  zmatval[s]=1; 
		  s++;
		}
	      
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,i,cqp->n)])&& i >= cqp->nb_int)
		{
		  /*y_ii >= 2u_ix_i - u_i^2*/ 
		  zmatind[s] =cqp->cont[cqp->ind_y_2 + i* (2*(cqp->n) -i +1)/2];
		  zmatval[s]=-1; 
		  s++;   
		  
		  /*y_ii >= 2l_ix_i - l_i^2*/ 
		  zmatind[s] =cqp->cont[cqp->ind_y_3 + i* (2*(cqp->n) -i +1)/2];
		  zmatval[s]=-1; 
		  s++;   
		  
		 
		}
	    }
	
	  if(j>i)
	    {	  
	     if (Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	       {
		  /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
		  zmatind[s] = cqp->cont[ i* (2*(cqp->n) -i +1)/2 + j-i];
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
		  zmatind[s] =cqp->cont[cqp->ind_y_2+ i* (2*(cqp->n) -i +1)/2 + j-i];
		  zmatval[s]=-1; 
		  s++;
		  
		  /*y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
		  zmatind[s] =cqp->cont[cqp->ind_y_3+ i* (2*(cqp->n) -i +1)/2 + j-i];
		  zmatval[s]=-1; 
		  s++;
		}
	    }
	}



   /* Second member of constrains */
   zrhs = alloc_vector_d (cqp->nb_row);
  i=0;
   
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=-(double)bab.u[k] *(double)bab.l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (Negative_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=-(double)bab.u[j] *(double)bab.l[k] ;
	  i++;
	}
   
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=(double)bab.u[k] *(double)bab.u[j] ;
	  i++;
	}
      
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (Positive_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  zrhs[i]=(double)bab.l[k] *(double)bab.l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    if ( (Negative_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]) && (j >= cqp->nb_int)) || (!Zero_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]) && (j < cqp->nb_int)))
       if (i<cqp->nb_row)
	 {
	   zrhs[i]=0;
	   i++;
	 }
  
  
  
  /* Sense of constraints */
  zsense    = alloc_string(cqp->nb_row);
   
  for(i=0;i<cqp->nb_row - cqp->nb_y_5;i++)
    zsense[i]= 'L';
   for(j=0;j<cqp->n;j++)
     if (Negative_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]) && (j >= cqp->nb_int))
       {
	 zsense[i]= 'L';
	 i++;
       }
     else
       if (!Zero_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]) && (j<cqp->nb_int))
	 {
	   zsense[i]= 'E';
	   i++;
	 }
   
  /***********************************************************************************/
  /*********************************** Bounds*****************************************/
  /***********************************************************************************/
 
  zlb= alloc_vector_d(cqp->nb_col);
  for(i=0;i<cqp->nb_col-1;i++)
    zlb[i]=(double)bab.l[i];
  zlb[cqp->nb_col-1] =  cqp->cons; 
  
  zub = alloc_vector_d(cqp->nb_col);
  for(i=0;i<cqp->nb_col-1;i++)
    zub[i]=(double)bab.u[i];
  zub[cqp->nb_col-1] = cqp->cons; 
  
  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
   
   
  /* Be sure that the new matrix is SDP, compute new lambda_min and add it to diagonal terms*/
  
 
  
  zobj = alloc_vector_d (cqp->nb_col);
  for(i=0;i<cqp->nb_col;i++)
    {
      if(i<cqp->n)
	zobj[i]= cqp->c[i];
      else
	zobj[i]=0;
    }

   
  i=cqp->ind_max_x;
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x ;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)])){
	if (j==k) 
	  {
	    zobj[i]= cqp->q[ij2k(j,k,cqp->n)] - cqp->new_q[ij2k(j,k,cqp->n)]/2 + cqp->new_lambda_min; 
	    i++;
	  }
	else
	  {
	    zobj[i]= 2*(cqp->q[ij2k(j,k,cqp->n)] - (cqp->new_q[ij2k(j,k,cqp->n)])/2);
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
   
 
   
  for(i=0;i<cqp->nb_col;i++)
    zctype[i]= 'C';
   
  strcpy (zprobname, "boxqp");
   
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



int eval_sbb_boxqp_cplex(struct miqcp_bab bab,  double ** sol_xy, double * inf_bound, int computeLocalSol)
{char     *probname = NULL;  
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
  double   sol_admissible;
  double fx;
 
  nodecount++;

  double * x_y = alloc_vector_d (cqp->nb_col);

  CPXENVptr    env = NULL;
  CPXLPptr      lp = NULL;
  int           status;
   
  int           cur_numrows, cur_numcols;

  /* Initialize the CPLEX environment */
  int k,l,temp;
     
  int *i =  alloc_vector(1);
  int *j = alloc_vector(1);

  
  /* update upper bounds*/
  refine_bounds(bab);
   
  /* load program into cplex*/
  env = CPXopenCPLEX (&status);
  char  errmsg[1024];
  if ( env == NULL ) 
    {
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return 0;
    }



  
  status = CPXsetintparam (env, 4012, 0);
  if ( status )
    {
      printf ("Failure to change param 4012, error %d.\n", status);
      return 0 ;
    }

  status = CPXsetdblparam (env, CPX_PARAM_EPGAP,REL_GAP);
  if ( status ) 
    {
      printf ("Failure to change relative gap, error %d.\n", status);
      return 0;
    }

   
  
  /* Set memory parameter*/
  status = CPXsetdblparam (env, CPX_PARAM_WORKMEM, MEMORY);
  if ( status )
    {
      printf ("Failure to change memory, error %d.\n", status);
      return 0;
    }

    
  
  /* Fill in the data for the problem.  */

  status = setqpproblemdata_sbb_boxqp_cplex (bab,&probname, &numcols, &numrows, &objsen, &obj, 
				     &rhs, &sense, &matbeg, &matcnt, &matind, &matval,
				       &qmatbeg, &qmatcnt, &qmatind, &qmatval,
				       &lb, &ub, &ctype);
   
  if ( status ) 
    {
      printf ( "Failed to build problem data arrays.\n");
      return 0;
    }

  /* Create the problem. */
  lp = CPXcreateprob (env, &status, probname);

  if ( lp == NULL ) 
    {
      printf ("Failed to create LP.\n");
      return 0;
    }
   
  /* Now copy the problem data into the lp */
  status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs, 
		      sense, matbeg, matcnt, matind, matval,
		      lb, ub, NULL);
   
  if ( status )
    {
      printf ( "Failed to copy problem data.\n");
      return 0;
    }
 status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
   if ( status ) {
      fprintf (stderr, "Failed to copy quadratic matrix.\n");
      return 0;
   }
  /* Optimize the problem and obtain solution. */
    

    status = CPXqpopt (env, lp);

    if ( status ) 
      {
	printf ( "Failed to optimize MIQP.\n");
	return 0;
      }

    
  solstat = CPXgetstat (env, lp);

  if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL && solstat !=CPXMIP_OPTIMAL_TOL )
      {
	printf ( "\n\nnot optimal, exit ...\n\nsolstat : %d \n\n",solstat);
	return 0;
      }

  /* Write the output to the screen. */
  status = CPXgetobjval (env, lp, &sol_admissible);
   
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      sol_admissible = MAX_SOL_BB;
      return 0;
    }
  

 
   cur_numrows = CPXgetnumrows (env, lp);
   cur_numcols = CPXgetnumcols (env, lp);

   
   status = CPXgetx (env, lp, x_y, 0, cur_numcols-1);
   if ( status ) 
     {
       printf ( "Failed to get optimal integer x.\n");
       return 0;
     }
   
   
   if (sol_admissible - ub_bb  > -EPS_BB )
     {
       printf("\nnode %ld\t current sol = %6lf > UB = %lf, branch pruned ",nodecount, sol_admissible,ub_bb);
       if (nodecount==1)
	 lb_bb = sol_admissible;
       return 0;
     }

   
   if (nodecount==1)
     lb_bb = sol_admissible;
   
   double gap;
   int is_int;
   
   /*Compute local sol*/
   if (computeLocalSol)
    {
      if (cqp->nb_int==0)
	{
	  compute_local_sol(&bab.l,&bab.u);
	  
	  fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	  fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	  fx = fx + cqp->cons;
	  if (fx - ub_bb < -EPS_BB)
	    {
	      ub_bb = fx;
	      for(k=0;k<cqp->n;k++)
		best_x[k]= x_y[k];
	    }
	  
	  if (cqp->local_sol_adm   - ub_bb < - EPS_BB)
	    {
	      ub_bb = cqp->local_sol_adm ;
	      for(k=0;k<cqp->n;k++)
		best_x[k] = cqp->local_sol[k];
	      
	      gap = (lb_bb - ub_bb);
	      gap = gap/ub_bb;
	      v_abs_ref(&gap);
	      
	      if (nodecount ==1)
		printf("\nnode %ld\t\tUB* (local solver): %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	      else
		printf("\nnode %ld\t\tUB* (local solver): %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	      
	    }
	  else
	    {
	      
	      gap = (lb_bb - ub_bb);
	      gap = gap/ub_bb;
	      v_abs_ref(&gap);
	      
	      if (nodecount ==1)
		printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	      else
		if (nodecount % PRINT_LEVEL_NODE == 0)
		  printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	      
	    }
	}
      else
	{
	  #ifdef HAVE_SCIP

	  status = compute_local_sol_scip_bb(&bab.l,&bab.u);
	  
	  fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	  fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	  fx = fx + cqp->cons;
	  if (fx - ub_bb < -EPS_BB)
	    {
	      ub_bb = fx;
	      for(k=0;k<cqp->n;k++)
		best_x[k]= x_y[k];
	    }
	  
	  if (status ==1)
		{
		  if (cqp->local_sol_adm   - ub_bb < - EPS_BB)
		    {
		      ub_bb = cqp->local_sol_adm ;
		      for(k=0;k<cqp->n;k++)
			best_x[k] = cqp->local_sol[k];
		      
		      gap = (lb_bb - ub_bb);
		      gap = gap/ub_bb;
		      v_abs_ref(&gap);
		      if (nodecount ==1)
			printf("\nnode %ld\t\tUB* (local solver): %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		      else
			printf("\nnode %ld\t\tUB* (local solver): %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
		      
		      
		    }
		  else
		    {
		      gap = (lb_bb - ub_bb);
		      gap = gap/ub_bb;
		      v_abs_ref(&gap);
		      
		      if (nodecount ==1)
			printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		      else
			if (nodecount % PRINT_LEVEL_NODE == 0)
			  printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
		      
		    }

		}
	  else
	    {  
		  gap = (lb_bb - ub_bb);
		  gap = gap/ub_bb;
		  v_abs_ref(&gap);
		  if (nodecount ==1)
		    printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		  else
		    if (nodecount % PRINT_LEVEL_NODE == 0)
		      printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
		  
		}
#endif
#ifndef HAVE_SCIP
	  is_int=1;
	  for(k=0;k<cqp->nb_int;k++)
	    if (!is_integer(x_y[k]))
	      {
		is_int=0;
		break;
	      }
	  
	  if (is_int==1)
	    {
	      fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	      fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	      fx = fx + cqp->cons;
	      if (fx - ub_bb < -EPS_BB)
		{
		  ub_bb = fx;
		  for(k=0;k<cqp->n;k++)
		    best_x[k]= x_y[k];
		  
		  
		  gap = (lb_bb - ub_bb);
		  gap = gap/ub_bb;
		  v_abs_ref(&gap);
		  
		  if (nodecount ==1)
		    printf("\nnode %ld\t\tUB*\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		  else
		    printf("\nnode %ld\t\tUB*\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
		}
	      else
		{
		  gap = (lb_bb - ub_bb);
		  gap = gap/ub_bb;
		  v_abs_ref(&gap);
		      
		  if (nodecount ==1)
		    printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		  else
		    if (nodecount % PRINT_LEVEL_NODE == 0)
		      printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
		  
		}
	      
	    }
	  else
	    {
	      gap = (lb_bb - ub_bb);
	      gap = gap/ub_bb;
	      v_abs_ref(&gap);
	      
	      if (nodecount ==1)
		printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	      else
		if (nodecount % PRINT_LEVEL_NODE == 0)
		  printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	      
	    }
#endif
	}
    }
   else
     {
       is_int=1;
       for(k=0;k<cqp->nb_int;k++)
	 if (!is_integer(x_y[k]))
	   {
	     is_int=0;
	     break;
	   }
       if (is_int==1)
	{
	  fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	  fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	  fx = fx + cqp->cons;
	  if (fx - ub_bb < -EPS_BB)
	    {
	      ub_bb = fx;
	      for(k=0;k<cqp->n;k++)
		best_x[k]= x_y[k];
	    
	  
	      gap = (lb_bb - ub_bb);
	      gap = gap/ub_bb;
	      v_abs_ref(&gap); 
	      if (nodecount ==1)
		printf("\nnode %ld\t\tUB*\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	      else
		printf("\nnode %ld\t\tUB*\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	    }
	  else
	    {
	      gap = (lb_bb - ub_bb);
	      gap = gap/ub_bb;
	      v_abs_ref(&gap); 
	      if (nodecount ==1)
		printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	      else
		if (nodecount % PRINT_LEVEL_NODE == 0)
		  printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	    }
	}
       else
	 {
	   gap = (lb_bb - ub_bb);
	   gap = gap/ub_bb;
	   v_abs_ref(&gap); 
	   if (nodecount ==1)
	     printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
	   else
	     if (nodecount % PRINT_LEVEL_NODE == 0)
	       printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
	}
     }
   
   
    /* clean list of leaves*/
  if (liste != NULL && computeLocalSol)
    clean_liste(liste, ub_bb);
  
  /*end Compute local sol*/
  
   
  
  /* Free up the problem as allocated by CPXcreateprob, if necessary */
  if ( lp != NULL )  
    {
      status = CPXfreeprob (env, &lp);
      if ( status )
	{
	  printf ( "CPXfreeprob failed, error code %d.\n", status);
	  return 0;
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
	  return 0;
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


  *inf_bound=sol_admissible;
  *sol_xy = x_y;

  return 1;

}

  
void sbb_boxqp_cplex( struct miqcp_bab bab)
{
 
  double   sol_admissible, fx;
  double * x_y = alloc_vector_d (cqp->nb_col);
  int k,m,p,ind;
  int *i =  alloc_vector(1);
  int *j = alloc_vector(1);
  
  struct miqcp_bab new_bab_1;
  struct miqcp_bab new_bab_2; 
  father=-1;
  
  long test_time;
  
  lb_bb = MIN_SOL_BB;
  
  /*last parameter: 1 means that it is necessary to compute a local sol at this node else 0 */
  int status = eval_sbb_boxqp_cplex(bab, &x_y,&sol_admissible,1);
if(status ==1)
  {
    insert_liste(liste,sol_admissible, x_y, bab, cqp->nb_col);
    lb_bb = liste->first->inf_bound_value;
  }
 else
   {
     if (status !=0)
       printf("\nproblem unfeasible\n");
     return;
   }
 
 int is_not_feasible;
  /*decides if a local sol has to be computed at this node*/
  int x_inf_new_u=0;
  double relative_gap;
  double absolute_gap;

  
 while(liste->size!=0)
    {

      test_time = time(NULL);
      if (test_time - start_time > TIME_LIMIT_BB)
	{
	  printf( "\n Time limit exeeded, time > %d \n",  TIME_LIMIT_BB);
	  break;
	}
      
      absolute_gap = (lb_bb -ub_bb);
      v_abs_ref(&absolute_gap);
      if (cqp->n == cqp->nb_int && absolute_gap < EPS_ABS_GAP)
	{
	  printf( "\nnode %ld\tBranch pruned, absolute_gap: %lf < %lf \n", nodecount,absolute_gap,EPS_ABS_GAP);
     	  break;
	  
	}
      
      relative_gap = absolute_gap/ub_bb;
      v_abs_ref(&relative_gap);
      if (relative_gap  < EPS_BB )
	{
	  printf( "\nnode %ld\tBranch pruned, relative_gap: %lf < %lf \n", nodecount,relative_gap,EPS_BB);
	  break;
	}
      
      lb_bb = liste->first->inf_bound_value;
      
      is_not_feasible = select_i_j_miqcrq(liste->first->sol_x, *liste->first->bab , i, j);
      if( is_not_feasible)
	{
	  if (cqp->nb_int > *i) 
	    {
	      new_bab_1 = update_variables_int_branch_1(*liste->first->bab,liste->first->sol_x,i,j);
	      
	      new_bab_2 = update_variables_int_branch_2(*liste->first->bab,liste->first->sol_x,i,j);
	      
	      suppress_first_liste(liste);
	      
	     if (test_bound_miqcp(new_bab_1,cqp->nb_col))
	       {
		 status = eval_sbb_boxqp_cplex(new_bab_1, &x_y, &sol_admissible,1);
		 if (status ==1)
		   {
		     if (ub_bb > sol_admissible)
		       {
			 absolute_gap = (sol_admissible -ub_bb);
			 relative_gap = absolute_gap/ub_bb;
			 v_abs_ref(&relative_gap);
			 if (relative_gap  > EPS_BB )
			   insert_liste(liste,sol_admissible, x_y, new_bab_1, cqp->nb_col);
		       }
		   }
	       }
	     
	     if (test_bound_miqcp(new_bab_2,cqp->nb_col))
	       {
		 status = eval_sbb_boxqp_cplex(new_bab_2, &x_y, &sol_admissible,1);
		    if (status ==1)
		      {
			if (ub_bb > sol_admissible)
			  {
			    absolute_gap = (sol_admissible -ub_bb);
			    relative_gap = absolute_gap/ub_bb;
			    v_abs_ref(&relative_gap);
			    if (relative_gap  > EPS_BB )
			      insert_liste(liste,sol_admissible, x_y, new_bab_2, cqp->nb_col);
			  }
		      }
	       }
	     
	     
	     if (liste->size!=0)
	       lb_bb=liste->first->inf_bound_value;
	    }
	  else
	    {
	      new_bab_1 = update_variables_cont_branch_1(*liste->first->bab, liste->first->sol_x,i,j);
	      new_bab_2 = update_variables_cont_branch_2(*liste->first->bab, liste->first->sol_x,i,j);
	      suppress_first_liste(liste);
	      
	      
	      if (x_y[*i] <  new_bab_1.u[*i])
		x_inf_new_u=1;
	      else
		x_inf_new_u=0;
	      
	      if (test_bound_miqcp(new_bab_1,cqp->nb_col))
		{
		  if (x_inf_new_u)
		    status = eval_sbb_boxqp_cplex(new_bab_1, &x_y, &sol_admissible,0);
		  else
		    status = eval_sbb_boxqp_cplex(new_bab_1, &x_y, &sol_admissible,1);
		  if (status ==1)
		    {
		      if (ub_bb > sol_admissible)
			{
			  absolute_gap = (sol_admissible -ub_bb);
			  relative_gap = absolute_gap/ub_bb;
			  v_abs_ref(&relative_gap);
			  if (relative_gap  > EPS_BB )
			    insert_liste(liste,sol_admissible, x_y, new_bab_1, cqp->nb_col);
			  
			}
		    }
		  
		}
	      
	      if (test_bound_miqcp(new_bab_2,cqp->nb_col))
		{
		  if  (x_inf_new_u)
		    status = eval_sbb_boxqp_cplex(new_bab_2, &x_y,& sol_admissible,1);
		  else
		    status = eval_sbb_boxqp_cplex(new_bab_2,  &x_y,& sol_admissible,0);
		  if (status ==1)
		    {
		      if (ub_bb > sol_admissible)
			{
			  absolute_gap = (sol_admissible -ub_bb);
			  relative_gap = absolute_gap/ub_bb;
			  v_abs_ref(&relative_gap);
			  if (relative_gap  > EPS_BB )
			    insert_liste(liste,sol_admissible, x_y, new_bab_2, cqp->nb_col);
			  
			}
		    }
		  
		}
	      if (liste->size!=0)
		lb_bb=liste->first->inf_bound_value;
	    }
	}
      else
	{
	  fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	  fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	  fx = fx + cqp->cons;
	  if ((ub_bb > fx) && (lb_bb < fx))
	    {
	      ub_bb =  fx;
	     for(k=0;k<cqp->n;k++)
	       best_x[k]= x_y[k];
	     
	     absolute_gap = (lb_bb -ub_bb);
	     v_abs_ref(&absolute_gap);
	     
	     relative_gap = absolute_gap/ub_bb;
	     v_abs_ref(&relative_gap);
	     
	     if (nodecount ==1)
	       printf("\nnode %ld\t\tUB* (bb leaf): %6lf\t\tLB: %6lf\t\trelative_gap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,relative_gap,time(NULL)-start_time);
	      else
		printf("\nnode %ld\t\tUB* (bb leaf): %6lf\t\tLB: %6lf\t\trelative_gap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,relative_gap,time(NULL)-start_time);
	     
	     
	     
	    }
	 
	 suppress_first_liste(liste);
	 if (liste != NULL)
	   clean_liste(liste,ub_bb);
       }
    }
 
 
 printf("\n\nEnd of branch-and-bound, all nodes evaluated\n");
  
}
  

 
     



  


