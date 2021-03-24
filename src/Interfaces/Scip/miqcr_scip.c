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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<unistd.h>
#include <time.h>
#include <string.h>


#include "utilities.h"
#include"quad_prog.h"
#include "liste.h"
#include"miqcr_scip.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

extern Q_MIQCP qqp;

extern long nodecount;

extern int TIME_LIMIT_BB;
extern int PRINT_LEVEL_SCIP;

extern double EPS_BB;
extern double EPS_BETA;
extern double EPS_LM;

extern double * best_x;
extern double ub_bb;

SCIP_VAR ** x;
SCIP_VAR ** obj;

static SCIP_RETCODE setupProblem_miqcr(SCIP* scip)
{
  char name[SCIP_MAXSTRLEN];
  SCIP_CONS* cons;
  SCIP_Real constante = -qqp->cons;
    
    
 
  int i,j,k,l,s,lim,lim_bis,lim_ter, s_bis;
  int cpt,cpt_bis;
  
  /* create empty problem */
  SCIP_CALL( SCIPcreateProbBasic(scip, "sbb_miqcr") );
 
  /* change to maximization if optimizing number of circles instead of rectangle area */
  SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   
  
  /*Allocation memoire du vecteur de variable*/
  
  SCIP_CALL( SCIPallocMemoryArray(scip, &x, qqp->nb_col) );

  /* variables x **/
  for( i = 0; i <qqp->nb_int; i++ )
    {
      
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, qqp->l[i], qqp->u[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }

  /* variables y **/
  for(k=0;k<qqp->nb_int;k++)
    for(l=0;l<qqp->ind_max_x;l++)
   	if (!Zero_BB(qqp->nb_var_y[ij2k(k,l,qqp->n)]))
	  {
	    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
	    SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, qqp->l[i], qqp->u[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
	    SCIP_CALL( SCIPaddVar(scip, x[i]) );
	    i++;
	  }
  
  

  /* variables z **/
 
  for(j=0;j<qqp->nb_int;j++)
    {
      lim=log_base_2(qqp->u[j])+1;
      for(l=0;l<qqp->nb_int;l++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(j,l,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, 0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
		SCIP_CALL( SCIPaddVar(scip, x[i]) );
		i++;
	      }
	}
   
    }
  /*variables t */
 
  for(i;i<qqp->nb_col;i++)
    {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, 0.0,1.0, 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }

  /*variable obj for the objective function*/
  SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
  SCIP_CALL( SCIPcreateVarBasic(scip, &obj[0], "objective",-SCIPinfinity(scip) ,SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
  SCIP_CALL( SCIPaddVar(scip, obj[0]) );
  
  
  /* linear constraint: sum_i Ax_i = b  */
  for (j=0;j<qqp->m;j++)
    {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name,0, NULL, NULL, qqp->b[j],qqp->b[j]) );
      for (i=0;i<qqp->n;i++)
	if (!Zero_BB(qqp->a[ij2k(j,i,qqp->n)])){
	  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qqp->a[ij2k(j,i,qqp->n)] ) );
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
    
  /* linear constraint: sum_i Dx_i <= e  */
  for (j=0;j<qqp->p;j++)
    {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),qqp->e[j]) );
      for (i=0;i<qqp->n;i++)
	if (!Zero_BB(qqp->d[ij2k(j,i,qqp->n)])){
	  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qqp->d[ij2k(j,i,qqp->n)] ) );
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }

    
  /*Quadratic constraints*/
  /* x^TAqx = bq  */
  for (k=0;k<qqp->mq;k++)
    {
      l = qqp->ind_max_x;
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, qqp->bq[k], qqp->bq[k]) );
      for( i = 0; i < qqp->n;i++ )
	{
	  if (!Zero_BB(qqp->aq[ijk2l(k,0,i+1,qqp->n+1,qqp->n+1)]))
	    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*qqp->aq[ijk2l(k,0,i+1,qqp->n+1,qqp->n+1)] ));
	  for( j = i; j < qqp->n; j++ )
	    if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	      {
		if (i==j && !Zero_BB(qqp->aq[ijk2l(k,i+1,i+1,qqp->n+1,qqp->n+1)]))
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], qqp->aq[ijk2l(k,i+1,i+1,qqp->n+1,qqp->n+1)] ));
		else
		  if (!Zero_BB(qqp->aq[ijk2l(k,i+1,j+1,qqp->n+1,qqp->n+1)]))
		    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 2*qqp->aq[ijk2l(k,i+1,j+1,qqp->n+1,qqp->n+1)] ));
		l++;
	      }
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
    }

  /* x^TDqx <= eq  */
  for (k=0;k<qqp->pq;k++)
    {
      l = qqp->ind_max_x;
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), qqp->eq[k]) );
      for( i = 0; i < qqp->n;i++ )
	{
	  if (!Zero_BB(qqp->dq[ijk2l(k,0,i+1,qqp->n+1,qqp->n+1)]))
	    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*qqp->dq[ijk2l(k,0,i+1,qqp->n+1,qqp->n+1)] ));
	  for( j = i; j < qqp->n; j++ )
	    if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	      {
		if (i==j && !Zero_BB(qqp->dq[ijk2l(k,i+1,i+1,qqp->n+1,qqp->n+1)]))
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], qqp->dq[ijk2l(k,i+1,i+1,qqp->n+1,qqp->n+1)] ));
		else
		  if (!Zero_BB(qqp->dq[ijk2l(k,i+1,j+1,qqp->n+1,qqp->n+1)]))
		    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 2*qqp->dq[ijk2l(k,i+1,j+1,qqp->n+1,qqp->n+1)] ));
		l++;
	      }
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
    }
  /* x_i = sum_k 2^k t_ik*/
  l=qqp->ind_max_x +qqp->nb_y+qqp->nb_z;
  for(i=0;i<qqp->ind_max_int;i++)
    {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL,0,0) );
      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 1 ));
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<lim;j++)
	{
	  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -pow(2,j)));
	  l++;
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }


   /*y_ij = sum_k 2^k z_ijk */
  l=qqp->ind_max_x +qqp->nb_y;
  s = qqp->ind_max_x;
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->nb_int;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    {
	      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL,0,0) );
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[s], 1 ));
	      s++;
	      for(k=0;k<lim;k++)
		{
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -pow(2,k)));
		  l++;
		}
	      SCIP_CALL( SCIPaddCons(scip, cons) );
	      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	    }
	}
    }

    
  l = qqp->ind_max_x;	  
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
  for(i=0;i<qqp->ind_max_x;i++)
    for(j=0;j<qqp->ind_max_x;j++)
      {
       if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	 {
	  
	   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), qqp->u[i] *qqp->u[j]));
	  if (i==j)
	    {
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*qqp->u[i]));
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	     
	    }
	  else
	    {
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], qqp->u[j]) );
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], qqp->u[i]) );
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	     
	    }
	  SCIP_CALL( SCIPaddCons(scip, cons) );
	  SCIP_CALL( SCIPreleaseCons(scip, &cons) ); 
	  l++;
	 }
       
      }
		  
  l = qqp->ind_max_x;	
  /* y_ii >= x_i*/   
  for(i=0;i<qqp->ind_max_x;i++)
    for(j=0;j<qqp->ind_max_x;j++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	{
	if(i==j && ((i < qqp->nb_int) ||  (i >= qqp->nb_int && qqp->l[i]>=1)))
	  {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
		l++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	  }
	else
	  {
	  if(i==j && (qqp->u[i]==1 && qqp->l[i]==0))
	    {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -0,0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], -1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
		l++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	    }
	  else
	    l++;
	    }
		 
	}

  int * sym  = alloc_vector(qqp->n*qqp->n);
  l = qqp->ind_max_x;
  for(i=0;i<qqp->ind_max_x;i++)
    {
      for(j=0;j<qqp->ind_max_x;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	  {
	    sym[i*qqp->n + j] = l;
	    l++;
	  }
	else
	  sym[i*qqp->n + j] = -1;
    }
  
  l = qqp->ind_max_x;
  /* y_ij == y_ji*/   
  for(i=0;i<qqp->ind_max_x;i++)
    {
      for(j=0;j<qqp->ind_max_x;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	  {
	    if (i<j)
	      {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0,0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons,  x[sym[j*qqp->n + i]], 1 ));
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
		l++;
	      }
	    else
	      l++;
	  }
    }
    
  
  
  /* z_ijk <= u_j_t_ik */
  cpt=0;
  
  l=qqp->ind_max_x +qqp->nb_y;

  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->nb_int;j++)
	{
	  s = qqp->ind_max_x+qqp->nb_y +qqp->nb_z;	  
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[cpt +s], -qqp->u[j]));
		l++;
		s++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	      }
	}
     
      cpt=cpt+lim;
    }
  
  /* z_ijk <= x_j*/		
  l=qqp->ind_max_x +qqp->nb_y;
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->nb_int;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], -1));
		l++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	      }
	}
    }
  
  /* z_ijk >= x_j + u_j_t_ik - u_j */
  cpt=0;
  l=qqp->ind_max_x +qqp->nb_y;
 
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->nb_int;j++)
	{
	  s = qqp->ind_max_x+qqp->nb_y +qqp->nb_z;
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),qqp->u[j] ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], 1));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[cpt + s], qqp->u[j]));
		l++;
		s++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	      }
	}
      cpt=cpt+ lim;
    }
    

  /* z_ijk >=0 */		
  l=qqp->ind_max_x +qqp->nb_y;
  for(i=0;i<qqp->nb_int;i++)
    {
      lim=log_base_2(qqp->u[i])+1;
      for(j=0;j<qqp->nb_int;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qqp->n)]))
	    for(k=0;k<lim;k++)
	      {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0 ));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1));
		l++;
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	      }
	}
    }

  
  /*objective function*/
  SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name,0, NULL, NULL, 0, NULL, NULL,NULL,constante,constante));
  SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, obj[0], -1.00) );
  for( i = 0; i < qqp->n; i++ )
    {
      SCIP_CALL(SCIPaddSquareCoefQuadratic(scip, cons, x[i],(qqp->new_q[ij2k(i,i,qqp->n)]- ((qqp->new_lambda_min)*2))/2));
      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], qqp->c[i] ) );
      for( j = i+1; j < qqp->n; j++ )
	{
	 SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j], qqp->new_q[ij2k(i,j,qqp->n)] )); 
      }
  }

  
  i=qqp->ind_max_x;
  for(j=0;j<qqp->ind_max_x;j++)
    for(k=0;k<qqp->ind_max_x ;k++)
      if (!Zero_BB(qqp->nb_var_y[ij2k(j,k,qqp->n)])){
	     if (j==k) 
	   {
		SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i],(qqp->q[ij2k(j,k,qqp->n)] - qqp->new_q[ij2k(j,k,qqp->n)]/2 + qqp->new_lambda_min) ) );
		i++;
	   }
	 else
	   {
		SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i],(qqp->q[ij2k(j,k,qqp->n)] - qqp->new_q[ij2k(j,k,qqp->n)]/2 + qqp->new_lambda_min) ));
		i++;
	   }
       
       

      }
  SCIP_CALL( SCIPaddCons(scip, cons) );
  SCIP_CALL( SCIPreleaseCons(scip, &cons) );

  
  return SCIP_OKAY;
}
	



int miqcr_scip()
{
		
		
   double   sol_admissible;
   double * x_y;
   double fx;
   
   
   SCIP* scip;
   int i,k,l;
   
   
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, qqp->nb_col) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
   
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   
  
   
   /* Fill in the data for the problem.  */
   
   
   SCIP_CALL( setupProblem_miqcr(scip));
   
   PRINT_LEVEL_SCIP =3;
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", EPS_BB) ); 
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", TIME_LIMIT_BB )); 
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel",5) );

   
   SCIP_CALL( SCIPsolve(scip) );
   
   
   SCIP_Real * xvals;
   SCIP_CALL( SCIPallocBufferArray(scip, &xvals, qqp->nb_col) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),qqp->nb_col,x, xvals) );
   
   SCIP_Real * objval; 
   SCIP_CALL( SCIPallocBufferArray(scip, &objval, 1) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),1,obj, objval) ); 

   nodecount =  SCIPgetNNodes(scip);
 
   // copy the solution in our structures
  
 
   x_y = alloc_vector_d(qqp->nb_col);
   for (i=0;i<qqp->nb_col;i++)
     x_y[i]=xvals[i];
   sol_admissible =objval[0];

   
   ub_bb  = sol_admissible ;
 
   
   int m;
   for(m=0;m<qqp->n;m++)
     best_x[m] = x_y[m];
   
   
   /* free memory arrays */
   for (i=0;i<qqp->nb_col;i++)
     SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
   
   SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
   
   SCIPfreeMemoryArray(scip, &obj);
   SCIPfreeMemoryArray(scip, &x);
   
   SCIP_CALL( SCIPfree(&scip) );
   
   return 1;
}
