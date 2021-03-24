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
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include<string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include"in_out.h"
#include"quad_prog.h"
#include "utilities.h"
#include"local_sol_ipopt.h"
#include "IpStdCInterface.h"

extern double MAX_SOL_BB;
extern int MAX_TIME_SOL_INIT;
extern int MAX_TIME_SOL_LOCAL;
extern double EPS_LOCAL_SOL;
extern int PRINT_LEVEL_IPOPT;

extern MIQCP qp;
extern C_MIQCP cqp;
extern int is_not_0_1;
extern int type;

/******************************************************************************/
/**************************LOCAL SOL BB **************************************/
/******************************************************************************/

Bool eval_f(Index n, Number* x, Bool new_x,
	    Number* obj_value, UserDataPtr user_data)
{
  assert (n == cqp->n);
  
  double fx=0;
  int i,j;
  for(i=0;i<cqp->n;i++)
    {
      fx=fx + cqp->c[i]*x[i];
      for(j=i;j<cqp->n;j++)
	if (i==j)
	  fx = fx + cqp->q[ij2k(i,j,cqp->n)]*x[i]*x[j];
	else
	  fx = fx + 2*cqp->q[ij2k(i,j,cqp->n)]*x[i]*x[j];
    }

   fx = fx + cqp->cons;
   *obj_value = fx;

 


  return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
		 Number* grad_f, UserDataPtr user_data)
{
  assert (n == cqp->n);
  int i,j;
    
  for(i=0;i<cqp->n;i++)
    {
      grad_f[i]=cqp->c[i];
      for(j=0;j<cqp->n;j++)
	grad_f[i]=grad_f[i] + 2*cqp->q[ij2k(i,j,cqp->n)]*x[j];


    }
  

 
  
  return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
	    Index m, Number* g, UserDataPtr user_data)
{
  assert(n == cqp->n);
  int tot_cont = cqp->m+cqp->p+cqp->mq+cqp->pq;
  assert(m == tot_cont);

  int k,i,j;

  
  for (k=0;k<cqp->m;k++)
    {
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	g[k] = g[k] + cqp->a[ij2k(k,i,cqp->n)]*x[i];
    }
    
  for (;k<cqp->m+cqp->p;k++)
    {
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	g[k] = g[k] + cqp->d[ij2k(k-cqp->m,i,cqp->n)]*x[i];
    }

  for (k;k<cqp->m+cqp->p+cqp->mq;k++)
    {
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	{
	  g[k] = g[k] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),0,i+1,cqp->n+1,cqp->n+1)]*x[i];
	  for (j=i;j<cqp->n;j++)
	    {
	      if (i==j)
		g[k] = g[k] + cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)]*x[i]*x[j];
	      else
		g[k] = g[k] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)]*x[i]*x[j];
	    }
	}
    }
 
    for (k;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
    {
      g[k] = 0;
      for (i=0;i<cqp->n;i++)
	{
	  g[k] = g[k] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),0,i+1,cqp->n+1,cqp->n+1)]*x[i];
	  for (j=i;j<cqp->n;j++)
	    {
	      if (i==j)
		g[k] = g[k] + cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)]*x[i]*x[j];
	      else
		g[k] = g[k] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)]*x[i]*x[j];
	    }
	}
    }
    
    return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
		Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *values,
		UserDataPtr user_data)
{
  int i,j,k;
  int l=0;
  if (values == NULL) {
    /* return the structure of the jacobian */

    /* this particular jacobian is dense */
    for (k=0;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
      {
	for (i=0;i<cqp->n;i++)
	  {
	    iRow[l] = k;
	    jCol[l] = i;
	    l++;
	  }
      }
    
    
  }
  else {
    /* return the values of the jacobian of the constraints */

    for (k=0;k<cqp->m;k++)
      {
	for (i=0;i<cqp->n;i++)
	  {
	    values[l] = cqp->a[ij2k(k,i,cqp->n)];
	    l++;
	  }
      }
    
    for (k;k<cqp->m+cqp->p;k++)
      {
	for (i=0;i<cqp->n;i++)
	  {
	    values[l] = cqp->d[ij2k(k-cqp->m,i,cqp->n)];
	    l++;
	  }
      }
    
    for (k;k<cqp->m+cqp->p+cqp->mq;k++)
      {
	for (i=0;i<cqp->n;i++)
	  {
	    values[l] = 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),0,i+1,cqp->n+1,cqp->n+1)];
	    for(j=0;j<cqp->n;j++)
	      {
		values[l] =values[l] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)]*x[j];
	      }
	    l++;
	  }

      }
    
    for (k;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
       {
	 for (i=0;i<cqp->n;i++)
	   {
	     values[l] = 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),0,i+1,cqp->n+1,cqp->n+1)];
	     for(j=0;j<cqp->n;j++)
	       {
		 values[l] =values[l] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)]*x[j];
		 
	       }
	    l++;
	   }
	  
       }

    
  
  }
  return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
	    Index m, Number *lambda, Bool new_lambda,
	    Index nele_hess, Index *iRow, Index *jCol,
	    Number *values, UserDataPtr user_data)
{
  Index idx = 0; /* nonzero element counter */
  Index row = 0; /* row counter for loop */
  Index col = 0; /* col counter for loop */
  if (values == NULL) {
    /* return the structure. This is a symmetric matrix, fill the lower left
     * triangle only. */

    /* the hessian for this problem is actually dense */
    idx=0;
    for (row = 0; row < cqp->n; row++) {
      for (col = row; col < cqp->n; col++) {
	iRow[idx] = row;

	jCol[idx] = col;

	idx++;
      }
    }
  
    
    assert(idx == nele_hess);
  }
  else {
    /* return the values. This is a symmetric matrix, fill the lower left
     * triangle only */

    /* fill the objective portion */
    int i,j;
    int l=0;
    for (i = 0; i < cqp->n; i++) {
      for (j = i; j < cqp->n; j++) {
	values[l] = 2*cqp->q[ij2k(i,j,cqp->n)];
	l++;
      }
    }


    
    /* add the portion for the quadratic constraints */
   
    int k=0;
    for(k=0;k<cqp->mq;k++)
      {
	l=0;
	for (i = 0; i < cqp->n; i++) {
	  for (j = i; j < cqp->n; j++) {
	    values[l] = values[l]+ lambda[k+cqp->m+cqp->p]*2*cqp->aq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)];
	    l++;
	  }
	}
      }

    for(k=0;k<cqp->pq;k++)
      {
	l=0;
	for (i = 0; i < cqp->n; i++) {
	  for (j = i; j < cqp->n; j++) {
	    values[l] = values[l]+ lambda[k+cqp->m+cqp->p+cqp->mq]*2*cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)];
	    l++;
	  }
	}
      }
   
    


}
  return TRUE;
}



int  compute_local_sol_ipopt(double * x_y, double ** l, double ** u) 
{
  Index n=cqp->n;                          /* number of variables */
  Index m=cqp->m+cqp->p+cqp->mq+cqp->pq;                          /* number of constraints  */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_g = NULL;               /* constraint multipliers     at the solution */
  Number* mult_x_L = NULL;             /* lower bound multipliers
					  at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
					  at the solution */
  Number obj;                          /* objective value */
  int i,j;                             /* generic counter */

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = cqp->n*m;
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  Index nele_hess = cqp->n*(cqp->n+1)/2;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
			    indices at 0 */

  /* set the number of variables and allocate space for the bounds */
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for (i=0; i<n; i++) {
    x_L[i] = l[0][i];
    x_U[i] = u[0][i];
  }

  /* set the number of constraints and allocate space for the bounds */
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  for (i=0;i<cqp->m;i++)
    {
      g_L[i] = cqp->b[i];
      g_U[i] = cqp->b[i];
    }
  
  for (;i<cqp->m+cqp->p;i++)
    {
      g_L[i] = -2e19;
      g_U[i] = cqp->e[i-cqp->m];
    }

  for (;i<cqp->m+cqp->p+cqp->mq;i++)
    {
      g_L[i] = cqp->bq[i - cqp->m- cqp->p];
      g_U[i] = cqp->bq[i- cqp->m-cqp->p];
    }
  
  for (;i<cqp->m+cqp->p+cqp->mq+cqp->pq;i++)
    {
      g_L[i] = -2e19;
      g_U[i] = cqp->eq[i - cqp->m- cqp->p-cqp->mq];
    }
 
  
  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
			   index_style, &eval_f, &eval_g, &eval_grad_f,
			   &eval_jac_g, &eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  
  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", EPS_LOCAL_SOL);
  /*AddIpoptStrOption(nlp, "mu_strategy", "adaptive");*/
  /*No output*/
  AddIpoptIntOption(nlp, "print_level",PRINT_LEVEL_IPOPT);
  /*Max cpu time*/
  AddIpoptNumOption(nlp, "max_cpu_time", MAX_TIME_SOL_LOCAL);
  /*AddIpoptStrOption(nlp, "output_file", "../data/local.sol");*/

  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  
  /* possibility to start from a point (for instance the best solution found)*/
   
  if (MAX_SOL_BB  == cqp->local_sol_adm)
    for (i=0; i<n; i++) 
      x[i] = x_y[i];
  else
    for (i=0; i<n; i++) 
      x[i] = random_double(x_L[i],x_U[i]);

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);
  mult_g = (Number*)malloc(sizeof(Number)*m);
  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, NULL);

  
  
  
   if (status == Solve_Succeeded) {
    /*copy local sol */
     
     for (i=0; i<n; i++) 
      cqp->local_sol[i] = x[i];
    cqp->local_sol_adm=obj;

    /* free allocated memory */
    FreeIpoptProblem(nlp);
  
    free(mult_x_L);
    free(mult_x_U);
    free(x_L);
    free(x_U);
    free(g_L);
    free(g_U);

    return 1;
   }
   else
     {
       
       /* free allocated memory */
       FreeIpoptProblem(nlp);
       
       free(mult_x_L);
       free(mult_x_U);
       free(x_L);
       free(x_U);
       free(g_L);
       free(g_U);
       return 0;
     }

 

}

/************************************************************************************************/
/**************  Compute local sol initialisation ***********************************************/
/************************************************************************************************/

Bool eval_f_init(Index n, Number* x, Bool new_x,
	    Number* obj_value, UserDataPtr user_data)
{
  assert (n == qp->n);
  
  double fx=0;
  int i,j;
  for(i=0;i<qp->n;i++)
    {
      fx=fx + qp->c[i]*x[i];
      for(j=i;j<qp->n;j++)
	if (i==j)
	  fx = fx + qp->q[ij2k(i,j,qp->n)]*x[i]*x[j];
	else
	  fx = fx + 2*qp->q[ij2k(i,j,qp->n)]*x[i]*x[j];
    }

   fx = fx + qp->cons;
  *obj_value = fx;


  return TRUE;
}

Bool eval_grad_f_init(Index n, Number* x, Bool new_x,
		 Number* grad_f, UserDataPtr user_data)
{
  assert (n == qp->n);
  int i,j;
    
  for(i=0;i<qp->n;i++)
    {
      grad_f[i]=qp->c[i];
      for(j=0;j<qp->n;j++)
	grad_f[i]=grad_f[i] + 2*qp->q[ij2k(i,j,qp->n)]*x[j];

    }

  
  
  return TRUE;
}

Bool eval_g_init(Index n, Number* x, Bool new_x,
	    Index m, Number* g, UserDataPtr user_data)
{
  assert(n == qp->n);
  int tot_cont = qp->m+qp->p+qp->mq+qp->pq;
  assert(m == tot_cont);

  int k,i,j;

  
  for (k=0;k<qp->m;k++)
    {
      g[k]=0;
      for (i=0;i<qp->n;i++)
	g[k] = g[k] + qp->a[ij2k(k,i,qp->n)]*x[i];
    }
    
  for (;k<qp->m+qp->p;k++)
    {
      g[k]=0;
      for (i=0;i<qp->n;i++)
	g[k] = g[k] + qp->d[ij2k(k-qp->m,i,qp->n)]*x[i];
    }

  for (k;k<qp->m+qp->p+qp->mq;k++)
    {
      g[k]=0;
      for (i=0;i<qp->n;i++)
	{
	  g[k] = g[k] + 2*qp->aq[ijk2l(k-(qp->m+qp->p),0,i+1,qp->n+1,qp->n+1)]*x[i];
	  for (j=i;j<qp->n;j++)
	    if (i==j)
		g[k] = g[k] + qp->aq[ijk2l(k-(qp->m+qp->p),i+1,j+1,qp->n+1,qp->n+1)]*x[i]*x[j];
	      else
		g[k] = g[k] + 2*qp->aq[ijk2l(k-(qp->m+qp->p),i+1,j+1,qp->n+1,qp->n+1)]*x[i]*x[j];
	}
    }
 
    for (k;k<qp->m+qp->p+qp->mq+qp->pq;k++)
    {
      g[k] = 0;
      for (i=0;i<qp->n;i++)
	{
	  g[k] = g[k] + 2*qp->dq[ijk2l(k-(qp->m+qp->p+qp->mq),0,i+1,qp->n+1,qp->n+1)]*x[i];
	  for (j=i;j<qp->n;j++)
	    {
	      if (i==j)
		g[k] = g[k] + qp->dq[ijk2l(k-(qp->m+qp->p+qp->mq),i+1,j+1,qp->n+1,qp->n+1)]*x[i]*x[j];
	      else
		g[k] = g[k] + 2*qp->dq[ijk2l(k-(qp->m+qp->p+qp->mq),i+1,j+1,qp->n+1,qp->n+1)]*x[i]*x[j]; 
	    }
	}
    }

    
    return TRUE;
}

Bool eval_jac_g_init(Index n, Number *x, Bool new_x,
		Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *values,
		UserDataPtr user_data)
{
  int i,j,k;
  int l=0;
  if (values == NULL) {
    /* return the structure of the jacobian */

    /* this particular jacobian is dense */
    for (k=0;k<qp->m+qp->p+qp->mq+qp->pq;k++)
      {
	for (i=0;i<qp->n;i++)
	  {
	    iRow[l] = k;
	    jCol[l] = i;
	    l++;
	  }
      }
  }
  else {
    /* return the values of the jacobian of the constraints */

    for (k=0;k<qp->m;k++)
      {
	for (i=0;i<qp->n;i++)
	  {
	    values[l] = qp->a[ij2k(k,i,qp->n)];
	    l++;
	  }
      }

    for (k;k<qp->m+qp->p;k++)
      {
	for (i=0;i<qp->n;i++)
	  {
	    values[l] = qp->d[ij2k(k-qp->m,i,qp->n)];
	    l++;
	  }
      }

    for (k;k<qp->m+qp->p+qp->mq;k++)
      {
	for (i=0;i<qp->n;i++)
	  {
	    values[l] = 2*qp->aq[ijk2l(k-(qp->m+qp->p),0,i+1,qp->n+1,qp->n+1)];
	    for(j=0;j<qp->n;j++)
	      values[l] =values[l] + 2*qp->aq[ijk2l(k-(qp->m+qp->p),i+1,j+1,qp->n+1,qp->n+1)]*x[j];
	    l++;
	  }
      }

    for (k;k<qp->m+qp->p+qp->mq+qp->pq;k++)
       {
	 for (i=0;i<qp->n;i++)
	   {
	     values[l] = 2*qp->dq[ijk2l(k-(qp->m+qp->p+qp->mq),0,i+1,qp->n+1,qp->n+1)];
	     for(j=0;j<qp->n;j++)
	       values[l] =values[l] + 2*qp->dq[ijk2l(k-(qp->m+qp->p+qp->mq),i+1,j+1,qp->n+1,qp->n+1)]*x[j];
	     l++;
	   }
       }

    
  }

  
  return TRUE;
}

Bool eval_h_init(Index n, Number *x, Bool new_x, Number obj_factor,
	    Index m, Number *lambda, Bool new_lambda,
	    Index nele_hess, Index *iRow, Index *jCol,
	    Number *values, UserDataPtr user_data)
{
  Index idx = 0; /* nonzero element counter */
  Index row = 0; /* row counter for loop */
  Index col = 0; /* col counter for loop */
  if (values == NULL) {
    /* return the structure. This is a symmetric matrix, fill the lower left
     * triangle only. */

    /* the hessian for this problem is actually dense */
    idx=0;
    for (row = 0; row < qp->n; row++) {
      for (col = row; col < qp->n; col++) {
	iRow[idx] = row;
	jCol[idx] = col;
	idx++;
      }
    }
  
    
    assert(idx == nele_hess);
  }
  else {
    /* return the values. This is a symmetric matrix, fill the lower left
     * triangle only */

    /* fill the objective portion */
    int i,j;
    int l=0;
    for (i = 0; i < qp->n; i++) {
      for (j = i; j < qp->n; j++) {
	values[l] = 2*qp->q[ij2k(i,j,qp->n)];
	l++;
      }
    }


    
    /* add the portion for the quadratic constraints */
   
    int k=0;
    for(k=0;k<qp->mq;k++)
      {
	l=0;
	for (i = 0; i < qp->n; i++) {
	  for (j = i; j < qp->n; j++) {
	    values[l] = values[l]+ lambda[k+qp->m+qp->p]*2*qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)];
	    l++;
	  }
	}
      }

    for(k=0;k<qp->pq;k++)
      {
	l=0;
	for (i = 0; i < qp->n; i++) {
	  for (j = i; j < qp->n; j++) {
	    values[l] = values[l]+ lambda[k+qp->m+qp->p+qp->mq]*2*qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)];
	    l++;
	  }
	}
      }
   
    
 
}
  return TRUE;
}



int  compute_local_sol_ipopt_init() 
{
  Index n=qp->n;                          /* number of variables */
  Index m=qp->m+qp->p+qp->mq+qp->pq;                          /* number of constraints  */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_g = NULL;               /* constraint multipliers     at the solution */
  Number* mult_x_L = NULL;             /* lower bound multipliers
					  at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
					  at the solution */
  Number obj;                          /* objective value */
  int i,j;                             /* generic counter */

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = qp->n*m;
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  Index nele_hess = qp->n*(qp->n+1)/2;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
			    indices at 0 */

  /* set the number of variables and allocate space for the bounds */
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for (i=0; i<n; i++) {
    x_L[i] = qp->l[i];
    x_U[i] = qp->u[i];
  }

  /* set the number of constraints and allocate space for the bounds */
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  for (i=0;i<qp->m;i++)
    {
      g_L[i] = qp->b[i];
      g_U[i] = qp->b[i];
    }
  
  for (;i<qp->m+qp->p;i++)
    {
      g_L[i] = -2e19;
      g_U[i] = qp->e[i-qp->m];
    }

  for (;i<qp->m+qp->p+qp->mq;i++)
    {
      g_L[i] = qp->bq[i - qp->m- qp->p];
      g_U[i] = qp->bq[i- qp->m-qp->p];
    }
  
  for (;i<qp->m+qp->p+qp->mq+qp->pq;i++)
    {
      g_L[i] = -2e19;
      g_U[i] = qp->eq[i - qp->m- qp->p-qp->mq];
    }
 
  
  
  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
			   index_style, &eval_f_init, &eval_g_init, &eval_grad_f_init,
			   &eval_jac_g_init, &eval_h_init);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  

  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
   AddIpoptNumOption(nlp, "tol", EPS_LOCAL_SOL);

   /*No output*/
  AddIpoptIntOption(nlp, "print_level",PRINT_LEVEL_IPOPT);

  /*Max cpu time*/
  AddIpoptNumOption(nlp, "max_cpu_time", MAX_TIME_SOL_INIT);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  
  
  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  
  /* possibility to start from a point (for instance the best solution found)*/
   
 
  for (i=0; i<n; i++) 
    x[i] = random_double(x_L[i],x_U[i]);
    
   /*else
     {
       for (i=0; i<n; i++) {
   	 x = x_L[i];
       }
       }*/

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);
  mult_g = (Number*)malloc(sizeof(Number)*m);
  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, NULL);

  
   if (status == Solve_Succeeded)
     {
       /*copy local sol */
       qp->local_sol = alloc_vector_d(qp->n);
	 for (i=0; i<n; i++) 
	   qp->local_sol[i] = x[i];
       qp->sol_adm=obj;
       /* free allocated memory */
       FreeIpoptProblem(nlp);
       free(x);
       free(mult_x_L);
       free(mult_x_U);
       free(x_L);
       free(x_U);
       free(g_L);
       free(g_U);
       return 1;
     }
   else
     {
       qp->local_sol = alloc_vector_d(qp->n);
       for (i=0; i<n; i++) 
	 qp->local_sol[i] = 0;
       qp->sol_adm=MAX_SOL_BB;

       /* free allocated memory */
       FreeIpoptProblem(nlp);
       free(x);
       free(mult_x_L);
       free(mult_x_U);
       free(x_L);
       free(x_U);
       free(g_L);
       free(g_U);
       return 0;
     }
 

  

}
