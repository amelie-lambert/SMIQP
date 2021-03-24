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
#include <time.h>
#include <stdlib.h>
#include <assert.h>


#include"quad_prog.h"
#include "utilities.h"
#include "liste.h"
#include "local_sol_ipopt.h"
#include "IpStdCInterface.h"
#include"sbb_miqcr_ipopt.h"

extern C_MIQCP cqp;
extern Liste liste;
extern SDP psdp;

extern int cont_lin; /* indicates if there is ony linear constraints */
extern int unconstrained;

extern long start_time;

extern double ub_bb;
extern  double lb_bb;
extern double * best_x;
extern long nodecount;

extern double MIN_SOL_BB;
extern double MAX_SOL_BB;
extern int TIME_LIMIT_BB;
extern int PRINT_LEVEL_IPOPT;
extern int PRINT_LEVEL_NODE;
extern double EPS_BB;
extern double EPS_BETA;
extern double EPS_LM;
extern double EPS_ABS_GAP;
extern double REL_GAP;
extern int is_not_0_1;


/* Function Implementations */
Bool eval_f_sbb_miqcr(Index n, Number* x, Bool new_x,
	    Number* obj_value, UserDataPtr user_data)
{
  assert (n == (cqp->nb_col-1));
  
  double fx=0;
  int i,j,k;

  /* quadratic terms **/
  for(i=0;i<cqp->n;i++)
    {
      fx=fx + cqp->c[i]*x[i];
      for(j=i;j<cqp->n;j++)
	if (i==j)
	  {
	    fx = fx + (cqp->new_q[ij2k(i,j,cqp->n)]- ((cqp->new_lambda_min)*2)) *x[i]*x[j]/2;
	  }
	else
	  {
	    fx = fx + cqp->new_q[ij2k(i,j,cqp->n)]*x[i]*x[j];
	  }
    }


  
  i=cqp->ind_max_x;
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x ;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)])){
	if (j==k) 
	  {
	    fx = fx + (cqp->q[ij2k(j,k,cqp->n)] - cqp->new_q[ij2k(j,k,cqp->n)]/2 + cqp->new_lambda_min)* x[i];
	    i++;
	  }
	else
	  {
	    fx= fx + 2*(cqp->q[ij2k(j,k,cqp->n)] - (cqp->new_q[ij2k(j,k,cqp->n)]/2))*x[i];
	    i++;
	  }
      }
  
   fx = fx + cqp->cons;
  *obj_value = fx;

  return TRUE;
}

Bool eval_grad_f_sbb_miqcr(Index n, Number* x, Bool new_x,
		 Number* grad_f, UserDataPtr user_data)
{
  assert (n == (cqp->nb_col-1));
  int i,j,k;
   
  for (k=0;k<cqp->nb_col-1;k++)
    grad_f[k]=0;
  
  for(i=0;i<cqp->n;i++)
    {
      grad_f[i]=cqp->c[i];
      for(j=0;j<cqp->n;j++)
	if (i==j)
	  {
	    grad_f[i]=grad_f[i] + (cqp->new_q[ij2k(i,j,cqp->n)]- ((cqp->new_lambda_min)*2))*x[j];
	  }
	  else
	    {
	      grad_f[i]=grad_f[i] + (cqp->new_q[ij2k(i,j,cqp->n)])*x[j];
	    }

    }
  
  i=cqp->ind_max_x;
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x ;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)])){
	if (j==k) 
	  {
	    grad_f[i] =  grad_f[i] + (cqp->q[ij2k(j,k,cqp->n)] - cqp->new_q[ij2k(j,k,cqp->n)]/2 + cqp->new_lambda_min); 
	    i++;
	  }
	else
	  {
	    grad_f[i] =  grad_f[i] + 2*(cqp->q[ij2k(j,k,cqp->n)] - (cqp->new_q[ij2k(j,k,cqp->n)]/2));
	    i++;
	  }

      }
 
  return TRUE;
}

Bool eval_g_sbb_miqcr(Index n, Number* x, Bool new_x,
	    Index m, Number* g, UserDataPtr user_data)
{
  assert(n == (cqp->nb_col-1));
  int tot_cont = cqp->nb_row;
  struct miqcp_bab* bab = user_data;
  assert(m == tot_cont);

  int k,i,j,l;

  
  for (k=0;k<m;k++)
    g[k]=0;

  k=0;
  /*Ax =b */
  for (k=0;k<cqp->m;k++)
    {
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	g[k] = g[k] + cqp->a[ij2k(k,i,cqp->n)]*x[i];
    }
  /*Dx <=e */ 
  for (k;k<cqp->m+cqp->p;k++)
    {
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	g[k] = g[k] + cqp->d[ij2k(k-cqp->m,i,cqp->n)]*x[i];
    }

  /* Aq Y = bq*/
 
  for (k;k<cqp->m+cqp->p+cqp->mq;k++)
    {
      l = cqp->ind_max_x;
      g[k]=0;
      for (i=0;i<cqp->n;i++)
	{
	  g[k] = g[k] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),0,i+1,cqp->n+1,cqp->n+1)]*x[i];
	  for (j=i;j<cqp->n;j++)
	    {
	     	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		  {
		    if (i==j)
		      g[k] = g[k] + cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)]*x[l];
		    else
		      g[k] = g[k] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)]*x[l];
		    l++;
		  }
	    }
	}
    }
  /* Dq Y <= eq*/

    for (k;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
    {
      l = cqp->ind_max_x;
      g[k] = 0;
      for (i=0;i<cqp->n;i++)
	{
	  g[k] = g[k] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),0,i+1,cqp->n+1,cqp->n+1)]*x[i];
	  for (j=i;j<cqp->n;j++)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		{
		  if (i==j)
		    g[k] = g[k] + cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)]*x[l];
		  else
		    g[k] = g[k] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)]*x[l];
		  l++;
		}
	    }
	}
    }

    
    l = cqp->ind_max_x;
    /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if(i==j)
	      g[k] = g[k]-((double)(*bab).l[i]+(double)(*bab).u[j])*x[i] + x[l];
	    else
	      g[k] = g[k] -(double)(*bab).l[i]*x[j] -(double)(*bab).u[j]*x[i] + x[l];
	    k++;
	    l++;
	  }

    /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
    
    l = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      l++;
	    else
	      {
		g[k] = g[k] -(double)(*bab).l[j]*x[i]-(double)(*bab).u[i]*x[j] + x[l];
		k++;
		l++;
	      }
	  }
	    
    /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
    
    l = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)])) 
	  {
	    if(i==j)
	      g[k] = g[k]+2*(double)(*bab).u[i]*x[i] - x[l];
	    else
	      g[k] = g[k] +(double)(*bab).u[i]*x[j] +  (double)(*bab).u[j]*x[i] - x[l];
	    k++;
	    l++;
	  }
    
    /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
    
    l = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if(i==j)
	      g[k] = g[k] + 2*(double)(*bab).l[i]*x[i] - x[l];
	    else
	      g[k] = g[k] + (double)(*bab).l[i]*x[j] + (double)(*bab).l[j]*x[i] - x[l];
	    k++;
	    l++;
	  }
	    
      
    /* y_ii >= x_i*/
    
    l = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for (j=i;j<cqp->n;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      {
		if((i < cqp->nb_int) ||  (i >= cqp->nb_int && cqp->l[i]>=1))
		  {
		    g[k] = g[k] + x[i] - x[l ];
		    k++;
		    l++;
		  }
		else
		  {
		    if(i >= cqp->nb_int && cqp->u[i]==1 && cqp->l[i]==0)
		      {
			g[k] = g[k] - x[i] + x[l ];
			k++;
			l++;
		      }
		    else
		      l++;
		  }
	      }
	    else
	      l++;
	  }
    
    
    return TRUE;
}



Bool eval_jac_g_sbb_miqcr(Index n, Number *x, Bool new_x,
		Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
  int i,j,k;
  int l=0;
  int ind_y;
  struct miqcp_bab* bab = user_data;
  if (values == NULL) {
    /* return the structure of the jacobian */

    /* this particular jacobian is dense */
    for (k=0;k<cqp->nb_row;k++)
      {
	for (i=0;i<(cqp->nb_col-1);i++)
	  {
	    iRow[l] = k;
	    jCol[l] = i;
	    l++;
	  }
      }
    
    
  }
  else {
    /* return the values of the jacobian of the constraints */
    int dim_values = cqp->nb_row*(cqp->nb_col-1);
    
    for (i=0;i<dim_values;i++)
      values[i]=0;


    /*Ax =b */
    for (k=0;k<cqp->m;k++)
      {
	for (i=0;i<(cqp->nb_col-1);i++)
	  {
	    if (i<cqp->ind_max_x)
	      values[l] =  cqp->a[ij2k(k,i,cqp->n)];
	    else
	      values[l]=0;
	    l++;
	  }
      }
    /*Dx <=e */ 
    for (k;k<cqp->m+cqp->p;k++)
      {
	for (i=0;i<(cqp->nb_col-1);i++)
	  {
	    if (i<cqp->ind_max_x)
	      values[l] = cqp->d[ij2k(k-cqp->m,i,cqp->n)];
	    else
	      values[l]=0;
	    l++;
	  }
      }
  
    /* Aq Y = bq*/
    for (k;k<cqp->m+cqp->p+cqp->mq;k++)
      {
	for (i=0;i<cqp->ind_max_x;i++)
	  {  
	    values[l] = 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),0,i+1,cqp->n+1,cqp->n+1)];
	    l++;
	  }
	for (i=0;i<cqp->ind_max_x;i++)
	  for(j=i;j<cqp->ind_max_x;j++)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		{
		  if (i==j)
		    values[l] =values[l] + cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)];
		  else
		    values[l] =values[l] + 2*cqp->aq[ijk2l(k-(cqp->m+cqp->p),i+1,j+1,cqp->n+1,cqp->n+1)];
		  l++;
		  }
	      
	    }
      }
    
   
    /* Dq Y <= eq*/
    for (k;k<cqp->m+cqp->p+cqp->mq+cqp->pq;k++)
      {
	for (i=0;i<cqp->ind_max_x;i++)
	  {
	    values[l] = 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),0,i+1,cqp->n+1,cqp->n+1)];
	    l++;
	  }
	
	for (i=0;i<cqp->ind_max_x;i++)
	  for(j=i;j<cqp->ind_max_x;j++)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
		{
		  if (i==j)
		    values[l] =values[l] + cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)];
		  else
		    values[l] =values[l] + 2*cqp->dq[ijk2l(k-(cqp->m+cqp->p+cqp->mq),i+1,j+1,cqp->n+1,cqp->n+1)];
		  l++;
		}
	      
	    }
      }
    
    
    l = (cqp->nb_col-1) * (cqp->m + cqp->p + cqp->mq+cqp->pq);
    /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
    ind_y = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      {
		values[l+i] =  - (double)(*bab).u[i] -(double)(*bab).l[i];
		values[l + ind_y] = 1;
		ind_y++;
	      }
	    else
	      {
		values[l+i] =  - (double)(*bab).u[j];
		values[l+j] =  -(double)(*bab).l[i];
		values[l + ind_y] = 1;
		ind_y++;
	      }
	    l = l + (cqp->nb_col-1);
	    
	  }
    
    /* y_ij <= l_jx_i + u_ix_j - l_ju_j*/
    ind_y = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      ind_y++;
	    else
	      {
		values[l+i] =  - (double)(*bab).l[j];
		values[l+j] =  -(double)(*bab).u[i];
		values[l + ind_y] = 1;
		ind_y++;
	      
		l = l + (cqp->nb_col-1);
	      }
	  }

    /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
    ind_y = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      {
		values[l+i] =  2*(double)(*bab).u[i];
		values[l + ind_y] = -1;
		ind_y++;
	      }
	    else
	      {
		values[l+i] =  (double)(*bab).u[j];
		values[l+j] =  (double)(*bab).u[i];
		values[l + ind_y] = -1;
		ind_y++;
	      }
	    l = l + (cqp->nb_col-1);
	  }
  
    
    /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
    ind_y = cqp->ind_max_x;
    for(i=0;i<cqp->ind_max_x;i++)
      for(j=i;j<cqp->ind_max_x;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (i==j)
	      {
		values[l+i] =  2*(double)(*bab).l[i];
		values[l + ind_y] = -1;
		ind_y++;
	      }
	    else
	      {
		values[l+i] =  (double)(*bab).l[j];
		values[l+j] =  (double)(*bab).l[i];
		values[l + ind_y] = -1;
		ind_y++;
	      }
	    l = l + (cqp->nb_col-1);
	  }
      
    /* y_ii >= x_i*/
      ind_y = cqp->ind_max_x;
      for(i=0;i<cqp->ind_max_x;i++)
	for(j=i;j<cqp->ind_max_x;j++)
	  if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	    {
	      if(i==j && ((i < cqp->nb_int) ||  (i >= cqp->nb_int && cqp->l[i]>=1)))
		{
		  values[l+i] = 1.0;
		  values[l+ind_y] = -1.0;
		  ind_y++;
		  l = l + (cqp->nb_col-1);
		}
	      else
		{
		  if(i==j && (i >= cqp->nb_int && cqp->u[i]==1 && cqp->l[i]==0))
		    {
		      values[l+i] = -1.0;
		      values[l+ind_y] = +1.0;
		      ind_y++;
		      l = l + (cqp->nb_col-1);
		    }
		  else
		    ind_y++;
		}
	     
	    }


   
  } 
  return TRUE;
}

Bool eval_h_sbb_miqcr(Index n, Number *x, Bool new_x, Number obj_factor,
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
	if (i==j)
	  values[l] = (cqp->new_q[ij2k(i,j,cqp->n)] - 2*cqp->new_lambda_min);
	else
	  values[l] = cqp->new_q[ij2k(i,j,cqp->n)]; 
	l++;
      }
    }

    
    


  }
  return TRUE;
}




/* Main Program */
int  eval_sbb_miqcr_ipopt(struct miqcp_bab bab,double ** sol_xy, double * sol_adm, int computeLocalSol) 
{
  Index n=(cqp->nb_col-1);                          /* number of variables */
  Index m=cqp->nb_row;                          /* number of constraints  */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x_y = NULL;                    /* starting point and solution vector */
  Number* mult_g = NULL;               /* constraint multipliers     at the solution */
  Number* mult_x_L = NULL;             /* lower bound multipliers
					  at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
					  at the solution */
  Number sol_admissible;                          /* objective value */
  int i,j,k,l;                             /* generic counter */

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = n*m;
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  Index nele_hess = cqp->n*(cqp->n+1)/2;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
			    indices at 0 */

  double fx;
  nodecount++;

  int feasible_for_const=1;

  // UserDataPtr user_data = (UserDataPtr)bab;
  
  /* update upper bounds*/
  int lowerbound=1;
  for(k=0;k<cqp->n;k++)
    if (cqp->l[k]!=0)
      {
  	lowerbound=0;
  	break;
      }
  
  if ((lowerbound==1) && (cont_lin ==0) &&  (unconstrained==0))
    refine_bounds(bab);



  /* set the number of variables and allocate space for the bounds */
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  k=0;
  l=0;

  for (i=0; i<cqp->ind_max_x; i++)
    {
      x_L[k] = (double)bab.l[l];
      x_U[k] = (double)bab.u[l];
      k++;
      l++;
    }
  
  for (i=0; i<cqp->ind_max_x; i++) 
    for (j=i; j<cqp->ind_max_x; j++)
      {
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    x_L[k] = (double)bab.l[l];
	    x_U[k] = (double)bab.u[l];
	    k++;
	 
	l++;
      }
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
 
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  g_L[i] = -2e19;
	  g_U[i] =-(double)bab.u[k] *(double)bab.l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j+1;k<cqp->ind_max_x;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	   g_L[i] = -2e19;
	  g_U[i] =-(double)bab.u[j] *(double)bab.l[k] ;
	  i++;
	}
   
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  g_L[i] = -2e19;
	  g_U[i] =(double)bab.u[k] *(double)bab.u[j] ;
	  i++;
	}
      
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{	  
	  g_L[i] = -2e19;
	  g_U[i] =(double)bab.l[k] *(double)bab.l[j] ;
	  i++;
	}
  
  for(j=0;j<cqp->ind_max_x;j++)
    if (!Zero_BB(cqp->nb_var_y[ij2k(j,j,cqp->n)]))
      {
	if( (j<cqp->nb_int && cqp->u[j]!=1 )  || (j >= cqp->nb_int && ((cqp->u[j]==1 && cqp->l[j]==0) ||  cqp->l[j]>=1)))
	  {	  
	    g_L[i] = -2e19;
	    g_U[i] =0;
	    i++;
	  }
	if(j < cqp->nb_int && cqp->u[j]==1 && cqp->l[j]==0)
	   {	  
	     g_L[i] =0;
	     g_U[i] =0;
	    i++;
	  }
      }


  
  
  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
			   index_style, &eval_f_sbb_miqcr, &eval_g_sbb_miqcr, &eval_grad_f_sbb_miqcr,
			   &eval_jac_g_sbb_miqcr, &eval_h_sbb_miqcr);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  
  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", REL_GAP);

  /*No output*/
  AddIpoptIntOption(nlp, "print_level",PRINT_LEVEL_IPOPT);
  
 

  /* Allocate space for the initial point and set the values */
  x_y = (Number*)malloc(sizeof(Number)*n);
  
  /* possibility to start from a point (for instance the best solution found)*/
   
  if (MAX_SOL_BB  != ub_bb)
    {
      for (i=0; i<cqp->n; i++) 
	x_y[i] = best_x[i];
      for (i;i<n;i++)
	x_y[i] = random_double(x_L[i],x_U[i]);
    }
  else
    for (i=0; i<n; i++) 
      x_y[i] = random_double(x_L[i],x_U[i]);

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);
  mult_g = (Number*)malloc(sizeof(Number)*m);
  /* solve the problem */
  status = IpoptSolve(nlp, x_y, NULL, &sol_admissible, mult_g, mult_x_L, mult_x_U, &bab);

   
  if (status != Solve_Succeeded) {
    
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
	status = compute_local_sol_ipopt(x_y,&bab.l,&bab.u);

      feasible_for_const= feasible_for_constraints(x_y);
  
      
      if (feasible_for_const)
	{
	  fx =sum_ij_qi_qj_x_ij(cqp->q, x_y, cqp->n);
	  fx =fx + sum_i_ci_x_i(cqp->c, x_y, cqp->n);
	  fx = fx + cqp->cons;
	  if (fx - ub_bb < -EPS_BB)
	    {
	      ub_bb = fx;
	      for(k=0;k<cqp->n;k++)
		best_x[k]= x_y[k];
	    }
	}

      
      if (status == 1)
	{
	  if (cqp->local_sol_adm  - ub_bb < -EPS_BB)
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
	  feasible_for_const= feasible_for_constraints(x_y);
	  if (feasible_for_const)
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
			printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,gap,time(NULL)-start_time);
		      else
			if (nodecount % PRINT_LEVEL_NODE == 0)
			  printf("\nnode %ld\t\tUB\t\t  : %6lf\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,gap,time(NULL)-start_time);
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
  
  *sol_adm=sol_admissible;
  *sol_xy = x_y;

 
  /* clean list of leaves*/

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


void sbb_miqcr_ipopt( struct miqcp_bab bab)
{

  double   sol_admissible,fx;
  double * x_y = alloc_vector_d (cqp->nb_col-1);

  int *i =  alloc_vector(1);
  int *j = alloc_vector(1);

  int k,m,p,ind;
  
  struct miqcp_bab new_bab_1; 
  struct miqcp_bab new_bab_2; 
  int test_time;
   int feasible_for_const;
  lb_bb = MIN_SOL_BB;
 
  

  /*last parameter: 1 means that it is necessary to compute a local sol at this node else 0 */
  int status = eval_sbb_miqcr_ipopt(bab,  &x_y,&sol_admissible,1);
  if(status ==1)
    {
      insert_liste(liste,sol_admissible, x_y, bab, (cqp->nb_col-1));
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
	     new_bab_1 = update_variables_int_branch_1(*liste->first->bab, liste->first->sol_x,i,j);
		     
	     new_bab_2 = update_variables_int_branch_2(*liste->first->bab,liste->first->sol_x,i,j);
		     
	     suppress_first_liste(liste);
	     
	     if (test_bound_miqcp(new_bab_1,(cqp->nb_col-1)))
	       {
		 status = eval_sbb_miqcr_ipopt(new_bab_1, &x_y, &sol_admissible,1);
	
		 if (status ==1)
		   {
		     if (ub_bb > sol_admissible)
		       {
			 absolute_gap = (sol_admissible -ub_bb);
			 relative_gap = absolute_gap/ub_bb;
			 v_abs_ref(&relative_gap);
			
			 if (relative_gap  > EPS_BB )
			   insert_liste(liste,sol_admissible, x_y, new_bab_1, (cqp->nb_col-1));

			
		       }
		   }
	       }
	     
	     if (test_bound_miqcp(new_bab_2,(cqp->nb_col-1)))
	       {
		 if (ub_bb > sol_admissible)
		   {
		     status = eval_sbb_miqcr_ipopt(new_bab_2, &x_y, &sol_admissible,1);
		     if (status ==1)
		       {
			 absolute_gap = (sol_admissible -ub_bb);
			 relative_gap = absolute_gap/ub_bb;
			 v_abs_ref(&relative_gap);
			 
			 if (relative_gap  > EPS_BB )
			   insert_liste(liste,sol_admissible, x_y, new_bab_2, (cqp->nb_col-1));
			
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
	     
	     if (test_bound_miqcp(new_bab_1,(cqp->nb_col-1)))
	       {
		 if (x_inf_new_u)
		   status=eval_sbb_miqcr_ipopt(new_bab_1,  &x_y, &sol_admissible,0);
		 else
		   status =eval_sbb_miqcr_ipopt(new_bab_1, &x_y, &sol_admissible,1);
		 if (status ==1)
		   {
		     if (ub_bb > sol_admissible)
		       {
			 absolute_gap = (sol_admissible -ub_bb);
			 relative_gap = absolute_gap/ub_bb;
			 v_abs_ref(&relative_gap);
			 if (relative_gap  > EPS_BB )
			   insert_liste(liste,sol_admissible, x_y, new_bab_1, (cqp->nb_col-1));
			 
		       }
		   }
	       }
	     
	     if (test_bound_miqcp(new_bab_2,(cqp->nb_col-1)))
	       {
		 if  (x_inf_new_u)
		   status = eval_sbb_miqcr_ipopt(new_bab_2,  &x_y,& sol_admissible,1);
		 else
		   status = eval_sbb_miqcr_ipopt(new_bab_2, &x_y,& sol_admissible,0);
		 if (status ==1)
		   {
		     if (ub_bb > sol_admissible)
		       {
			 absolute_gap = (sol_admissible -ub_bb);
			 relative_gap = absolute_gap/ub_bb;
			 v_abs_ref(&relative_gap);
			 if (relative_gap  > EPS_BB )
			   insert_liste(liste,sol_admissible, x_y, new_bab_2, (cqp->nb_col-1));
			
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
	 fx =fx + cqp->cons;
	 feasible_for_const= feasible_for_constraints(x_y);
	 if ((ub_bb > fx) && (lb_bb < fx) && feasible_for_const)
	   {
	     ub_bb =  fx;
	     for(k=0;k<cqp->n;k++)
	       best_x[k]= x_y[k];
	      
	     absolute_gap = (lb_bb -ub_bb);
	     v_abs_ref(&absolute_gap);
	     
	     relative_gap = absolute_gap/ub_bb;
	     v_abs_ref(&relative_gap);
	     
	     if (nodecount ==1)
	       printf("\nnode %ld\t\tUB* (bb leaf)\t  : %6lf\t\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, sol_admissible,relative_gap,time(NULL)-start_time);
	      else
		printf("\nnode %ld\t\tUB* (bb leaf)\t  : %6lf\t\t\tLB: %6lf\t\tgap: %6lf\t\ttime:%ld",nodecount,ub_bb, lb_bb,relative_gap,time(NULL)-start_time);
	     
	     
	    
	   }
	 
	 suppress_first_liste(liste);
	 if (liste != NULL)
	   clean_liste(liste,ub_bb);
       }
    }
  printf("\n\nEnd of branch-and-bound, all nodes evaluated\n");

}





