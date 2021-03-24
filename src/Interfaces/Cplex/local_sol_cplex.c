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
#include"quad_prog.h"
#include "utilities.h"
#include"local_sol_cplex.h"

extern MIQCP qp;
extern C_MIQCP cqp;
extern int is_not_0_1;
extern int type;
extern double EPS_LOCAL_SOL;
extern int MAX_TIME_SOL_INIT;
extern int MAX_TIME_SOL_LOCAL;


/* parameters for CPLEX */
extern int MEMORY;
extern  double REL_GAP;
extern double * best_x;
extern double ub_bb;

int setqpproblemdata_cplex_local (double ** l, double ** u, char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p,double **lb_p, double **ub_p)

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
 
   
  
   
  int i,j,k,s;
   
  zprobname = alloc_string(16); 
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
 
  zmatbeg   = alloc_vector(cqp->n);   
  zmatcnt   = alloc_vector(cqp->n); 
    
  
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i< cqp->n;i++)
    zmatbeg[i]= zmatbeg[i-1] + cqp->m + cqp->p;
   
  
  int numnz;
  
  numnz=zmatbeg[cqp->n-1] + cqp->m + cqp->p;
  
    
  for(i=0;i<cqp->n-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  zmatcnt[cqp->n - 1] = numnz - zmatbeg[cqp->n-1];
 
 /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
  /* Contraints over x (order x_1 x_2 ... x_n ) */
 
  for(i=0;i<cqp->n;i++)
    {
      /*Ax = b*/
      for(j= 0;j<cqp->m;j++)
	{ 
	  zmatval[s]=cqp->a[ij2k(j,i,cqp->n)];
	  zmatind[s]=j;
	  s++;
	}
       
      /* Dx =e*/
      for(j=0;j<cqp->p;j++)
	{
	  zmatval[s]=cqp->d[ij2k(j,i,cqp->n)];
	  zmatind[s]=cqp->m + j;
	  s++;
	}
    }
  
  /* Second member of constrains */ 
  zrhs = alloc_vector_d (cqp->m + cqp->p);
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
    
  /* Sense of constraints */
  zsense = alloc_string(cqp->m + cqp->p);
   
  for(i=0;i<cqp->m;i++)
    zsense[i]= 'E';
  for(;i<cqp->m + cqp->p;i++)
    zsense[i]= 'L';

  
  
  /***********************************************************************************/
  /****************************** Bounds on variables ********************************/
  /***********************************************************************************/
   
  zlb= alloc_vector_d(cqp->n);
  for(i=0;i<cqp->n;i++)
    zlb[i]=l[0][i];
   
  zub = alloc_vector_d(cqp->n);
   
  for(i=0;i<cqp->n ;i++)
    zub[i]=u[0][i];
 

  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
 
  zobj = alloc_vector_d (cqp->n);
  for(i=0;i<cqp->n;i++)
    {
      if(i<cqp->n)
	    zobj[i]= cqp->c[i];
      else
	zobj[i]=0;
    }
 

  

  zqmatbeg  = alloc_vector(cqp->n);
  zqmatbeg[0]=0;
  for (i=1;i<cqp->n;i++)
    zqmatbeg[i]=zqmatbeg[i-1]+ cqp->n;
   
   
  zqmatcnt   = alloc_vector(cqp->n);
  for (i=0;i<cqp->n;i++)
    zqmatcnt[i] = cqp->n;
   

  zqmatind   = alloc_vector(cqp->n*cqp->n);
  zqmatval   = alloc_vector_d(cqp->n*cqp->n);
   
  k=0;
  for(i=0;i<cqp->n;i++)
    for(j=0;j<cqp->n;j++)
      { 
	zqmatval[k] =  2*cqp->q[ij2k(i,j,cqp->n)]; 
	zqmatind[k] = j;
	k++;
      }

  
  strcpy (zprobname, "local_cplex");

  
   
  *numcols_p   = cqp->n;
  *numrows_p   = cqp->m +cqp->p;
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
  
  


   
  return (status);
} /* END setproblemdata */

void compute_local_sol(double ** l, double ** u)
{
/* Declare pointers for the variables and arrays that will contain
      the data which define the LP problem.  The setproblemdata() routine
      allocates space for the problem data.  */

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
   

   
   int           cur_numrows, cur_numcols,i;
   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat;
   double   objval;
   cqp->local_sol = alloc_vector_d (cqp->n);
   
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
      return;
   }

  

   /* Fill in the data for the problem.  */

    status = setqpproblemdata_cplex_local (l,u, &probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval, &qmatbeg, &qmatcnt, &qmatind, &qmatval,&lb, &ub);
  
    if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return;
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
      return;
   }

   /* Now copy the LP part of the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      return;
   }

   status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
   if ( status ) {
      fprintf (stderr, "Failed to copy quadratic matrix.\n");
      return;
   }


   /* When a non-convex objective function is present, CPLEX will
      return error CPXERR_Q_NOT_POS_DEF unless the parameter
      CPXPARAM_SolutionTarget is set to accept first-order optimal
      solutions.  */
   status = CPXsetintparam (env, CPXPARAM_OptimalityTarget, CPX_OPTIMALITYTARGET_FIRSTORDER);
   if ( status ) return;

  
     status = CPXsetdblparam (env, CPX_PARAM_TILIM,MAX_TIME_SOL_LOCAL);
  if ( status ) 
    {
      printf ("Failure to change time limit, error %d.\n", status);
      return;
    }

  /* Optimize the problem and obtain solution. */

   status = CPXqpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize QP.\n");
      return;
   }
  
   solstat = CPXgetstat (env, lp);
   
   if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL && solstat!=CPX_STAT_NUM_BEST && solstat!=CPX_STAT_FIRSTORDER)
     {
       printf ( "\n\nCompute local sol failed, exit ...\n\n solstat : %d \n\n",solstat);
       return;
     }
   
  
  status = CPXgetobjval (env, lp, &cqp->local_sol_adm);
  
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return;
    }

  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  /*status = CPXgetx (env, lp, bab.x_local, 0, cur_numcols-1);*/
  status = CPXgetx (env, lp, cqp->local_sol, 0, cur_numcols-1);
 
   
  if ( status ) 
    {
      printf ( "Failed to get optimal x.\n");
      return;
    }
   
 
  

  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);
  
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
  
  
  

}



/**********************************************************************************************/
/***********  Compute local sol for presolve  *************************************************/
/**********************************************************************************************/



int setqpproblemdata_cplex_local_init ( char **probname_p, int *numcols_p, int *numrows_p, int *objsen_p, double **obj_p, double **rhs_p, char **sense_p, int **matbeg_p, int **matcnt_p, int **matind_p, double **matval_p, int **qmatbeg_p, int **qmatcnt_p, int **qmatind_p, double **qmatval_p, double **lb_p, double **ub_p)

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
 
   
  
   
  int i,j,k,s;
   
  zprobname = alloc_string(16); 
  /***********************************************************************************/
  /*********************************** Variables *************************************/
  /***********************************************************************************/
 
  zmatbeg   = alloc_vector(qp->n);   
  zmatcnt   = alloc_vector(qp->n); 
    
  
  /* x */
  zmatbeg[0]= 0;
  for(i=1;i< qp->n;i++)
    zmatbeg[i]= zmatbeg[i-1] + qp->m + qp->p;
   
  
  int numnz;
  
  numnz=zmatbeg[qp->n-1] + qp->m + qp->p;
  
    
  for(i=0;i<qp->n-1;i++)
    zmatcnt[i] = zmatbeg[i+1]- zmatbeg[i];
  zmatcnt[qp->n - 1] = numnz - zmatbeg[qp->n-1];
 
 /***********************************************************************************/
  /*********************************** Constraints ***********************************/
  /***********************************************************************************/
  zmatind   = alloc_vector(numnz);    
  zmatval   = alloc_vector_d(numnz);

  s=0;
   
  /* Contraints over x (order x_1 x_2 ... x_n ) */
 
  for(i=0;i<qp->n;i++)
    {
      /*Ax = b*/
      for(j= 0;j<qp->m;j++)
	{ 
	  zmatval[s]=qp->a[ij2k(j,i,qp->n)];
	  zmatind[s]=j;
	  s++;
	}
       
      /* Dx =e*/
      for(j=0;j<qp->p;j++)
	{
	  zmatval[s]=qp->d[ij2k(j,i,qp->n)];
	  zmatind[s]=qp->m + j;
	  s++;
	}
    }
  
  /* Second member of constrains */ 
  zrhs = alloc_vector_d (qp->m + qp->p);
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
    
  /* Sense of constraints */
  zsense = alloc_string(qp->m + qp->p);
   
  for(i=0;i<qp->m;i++)
    zsense[i]= 'E';
  for(;i<qp->m + qp->p;i++)
    zsense[i]= 'L';

  
  
  /***********************************************************************************/
  /****************************** Bounds on variables ********************************/
  /***********************************************************************************/
   
  zlb= alloc_vector_d(qp->n);
  for(i=0;i<qp->n;i++)
    zlb[i]=qp->l[i];
   
  zub = alloc_vector_d(qp->n);
   
  for(i=0;i<qp->n ;i++)
    zub[i]=qp->u[i];
 

  /**************************************************************************************/
  /******************************* Objective function ***********************************/
  /**************************************************************************************/
 
  zobj = alloc_vector_d (qp->n);
  for(i=0;i<qp->n;i++)
    {
      if(i<qp->n)
	    zobj[i]= qp->c[i];
      else
	zobj[i]=0;
    }
 

  

  zqmatbeg  = alloc_vector(qp->n);
  zqmatbeg[0]=0;
  for (i=1;i<qp->n;i++)
    zqmatbeg[i]=zqmatbeg[i-1]+ qp->n;
   
   
  zqmatcnt   = alloc_vector(qp->n);
  for (i=0;i<qp->n;i++)
    zqmatcnt[i] = qp->n;
   

  zqmatind   = alloc_vector(qp->n*qp->n);
  zqmatval   = alloc_vector_d(qp->n*qp->n);
   
  k=0;
  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      { 
	zqmatval[k] =  2*qp->q[ij2k(i,j,qp->n)]; 
	zqmatind[k] = j;
	k++;
      }

  
  strcpy (zprobname, "local_cplex");

  
   
  *numcols_p   = qp->n;
  *numrows_p   = qp->m +qp->p;
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
  
  


   
  return (status);
} /* END setproblemdata */

void compute_local_sol_init()
{
/* Declare pointers for the variables and arrays that will contain
      the data which define the LP problem.  The setproblemdata() routine
      allocates space for the problem data.  */

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
   int      *qmatbeg = NULL;
   int      *qmatcnt = NULL;
   int      *qmatind = NULL;
   double   *qmatval = NULL;

   
   int           cur_numrows, cur_numcols,i;
   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat;
   double   objval;
   qp->local_sol = alloc_vector_d (qp->n);
   
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
      return;
   }

  
   

   /* Fill in the data for the problem.  */

   status = setqpproblemdata_cplex_local_init (&probname, &numcols, &numrows, &objsen, &obj, &rhs, &sense, &matbeg, &matcnt, &matind, &matval,&qmatbeg, &qmatcnt, &qmatind, &qmatval, &lb, &ub);
  
    if ( status ) {
      fprintf (stderr, "Failed to build problem data arrays.\n");
      return;
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
      return;
   }

   /* Now copy the LP part of the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      return;
   }

   status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
   if ( status ) {
      fprintf (stderr, "Failed to copy quadratic matrix.\n");
      return;
   }


   /* When a non-convex objective function is present, CPLEX will
      return error CPXERR_Q_NOT_POS_DEF unless the parameter
      CPXPARAM_SolutionTarget is set to accept first-order optimal
      solutions.  */
   status = CPXsetintparam (env, CPXPARAM_OptimalityTarget, CPX_OPTIMALITYTARGET_FIRSTORDER);
   if ( status ) return;

   /* Optimize the problem and obtain solution. */

  
     status = CPXsetdblparam (env, CPX_PARAM_TILIM,MAX_TIME_SOL_INIT);
  if ( status ) 
    {
      printf ("Failure to change time limit, error %d.\n", status);
      return;
    }
   status = CPXqpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize QP.\n");
      return;
   }

   solstat = CPXgetstat (env, lp);
   
   if (solstat !=   CPX_STAT_OPTIMAL && solstat !=   CPXMIP_OPTIMAL && solstat!=CPX_STAT_NUM_BEST && solstat!=CPX_STAT_FIRSTORDER)
     {
       printf ( "\n\nCompute local sol failed, exit ...\n\n solstat : %d \n\n",solstat);
       return;
     }
   
  
  status = CPXgetobjval (env, lp, &qp->sol_adm);
  
  if ( status ) 
    {
      printf ("No MIP objective value available.  Exiting...\n");
      return;
    }

  
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);

  /*status = CPXgetx (env, lp, bab.x_local, 0, cur_numcols-1);*/
  status = CPXgetx (env, lp, qp->local_sol, 0, cur_numcols-1);
  
  /*for (i=0;i<qp->n;i++)
    x_local[0][i] = qp->local_sol[i];*/
   
  if ( status ) 
    {
      printf ( "Failed to get optimal x.\n");
      return;
    }
   
 
  cur_numrows = CPXgetnumrows (env, lp);
  cur_numcols = CPXgetnumcols (env, lp);
  
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
   free_and_null ((char **) &qmatbeg);
   free_and_null ((char **) &qmatcnt);
   free_and_null ((char **) &qmatind);
   free_and_null ((char **) &qmatval);
  
  

}





