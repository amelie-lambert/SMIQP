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
#include "local_sol_scip.h"
#include "liste.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "sbb_miqcr_scip.h"

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
extern int father;

extern double MIN_SOL_BB;
extern double MAX_SOL_BB;
extern int TIME_LIMIT_BB;
extern int PRINT_LEVEL_SCIP;
extern int PRINT_LEVEL_NODE;
extern double EPS_BB;
extern double EPS_BETA;
extern double EPS_LM;
extern double EPS_ABS_GAP;



SCIP_VAR ** x;
SCIP_VAR ** obj;


static SCIP_RETCODE setupProblem_sbb_miqcr(SCIP* scip,struct miqcp_bab bab)
{
  char name[SCIP_MAXSTRLEN];
  SCIP_CONS* cons;
  SCIP_Real constante = -cqp->cons;
    
    
 
  int i, j, k, l;
 
  /* create empty problem */
  SCIP_CALL( SCIPcreateProbBasic(scip, "sbb_miqcr") );
 
  /* change to maximization if optimizing number of circles instead of rectangle area */
  SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   
  
  /*Allocation memoire du vecteur de variable*/
  
  SCIP_CALL( SCIPallocMemoryArray(scip, &x, cqp->nb_col -1) );
    
  
    
  for( i=0; i <cqp->nb_col -1; i++ )
    {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, bab.l[i], bab.u[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }
    
  /*variable obj for the objective function*/
  SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
  SCIP_CALL( SCIPcreateVarBasic(scip, &obj[0], "objective",-SCIPinfinity(scip) ,SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
  SCIP_CALL( SCIPaddVar(scip, obj[0]) );
  
  
  /* linear constraint: sum_i Ax_i = b  */
  for (j=0;j<cqp->m;j++)
    {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name,0, NULL, NULL, cqp->b[j],cqp->b[j]) );
      for (i=0;i<cqp->n;i++)
	if (!Zero_BB(cqp->a[ij2k(j,i,cqp->n)])){
	  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],cqp->a[ij2k(j,i,cqp->n)] ) );
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
    
  /* linear constraint: sum_i Dx_i <= e  */
  for (j=0;j<cqp->p;j++)
    {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),cqp->e[j]) );
      for (i=0;i<cqp->n;i++)
	if (!Zero_BB(cqp->d[ij2k(j,i,cqp->n)])){
	  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],cqp->d[ij2k(j,i,cqp->n)] ) );
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }

    
  /*Quadratic constraints*/
  /* x^TAqx = bq  */
  for (k=0;k<cqp->mq;k++)
    {
      l = cqp->ind_max_x;
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, cqp->bq[k], cqp->bq[k]) );
      for( i = 0; i < cqp->n;i++ )
	{
	  if (!Zero_BB(cqp->aq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]))
	    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*cqp->aq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)] ));
	  for( j = i; j < cqp->n; j++ )
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	      {
		if (i==j && !Zero_BB(cqp->aq[ijk2l(k,i+1,i+1,cqp->n+1,cqp->n+1)]))
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], cqp->aq[ijk2l(k,i+1,i+1,cqp->n+1,cqp->n+1)] ));
		else
		  if (!Zero_BB(cqp->aq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]))
		    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 2*cqp->aq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)] ));
		l++;
	      }
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
    }

  /* x^TDqx <= eq  */
  for (k=0;k<cqp->pq;k++)
    {
      l = cqp->ind_max_x;
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), cqp->eq[k]) );
      for( i = 0; i < cqp->n;i++ )
	{
	  if (!Zero_BB(cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]))
	    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)] ));
	  for( j = i; j < cqp->n; j++ )
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	      {
		if (i==j && !Zero_BB(cqp->dq[ijk2l(k,i+1,i+1,cqp->n+1,cqp->n+1)]))
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], cqp->dq[ijk2l(k,i+1,i+1,cqp->n+1,cqp->n+1)] ));
		else
		  if (!Zero_BB(cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]))
		    SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 2*cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)] ));
		l++;
	      }
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
    }
  
  l = cqp->ind_max_x;
  /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      {
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),-(double)bab.u[j] *(double)bab.l[i]));
	    if (i==j)
	      {
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], -(double)bab.l[i]-(double)bab.u[j]));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
	      }
	    else
	      {
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], -(double)bab.u[j]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], -(double)bab.l[i]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
	      }
	      SCIP_CALL( SCIPaddCons(scip, cons) );
	      SCIP_CALL( SCIPreleaseCons(scip, &cons) ); 
	      l++;
	  }
      }
  
		  
  l = cqp->ind_max_x;	    
  /* y_ij <= l_jx_i + u_ix_j - l_ju_i*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	{
	if(i==j)
	   l++;
	else
	  {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), -(double)bab.u[i] *(double)bab.l[j]));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i],- (double)bab.l[j]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], -(double)bab.u[i]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
		SCIP_CALL( SCIPaddCons(scip, cons) );
		SCIP_CALL( SCIPreleaseCons(scip, &cons) );
		l++;
	   }
	}
		

  l = cqp->ind_max_x;	  
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      {
       if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	 {
	  SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), (double)bab.u[i] *(double)bab.u[j]));
	  if (i==j)
	    {
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*(double)bab.u[i]));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	    }
	  else{
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], (double)bab.u[j]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], (double)bab.u[i]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	      }
	  SCIP_CALL( SCIPaddCons(scip, cons) );
	  SCIP_CALL( SCIPreleaseCons(scip, &cons) ); 
	  l++;
	 }
     }
		  
  l = cqp->ind_max_x;			
  /* y_ij >= l_ix_j + l_jx_i - l_il_j*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      {
       if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	 {
	  SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),(double)bab.l[i] *(double)bab.l[j] ));
	  if (i==j)
	    {
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 2*(double)bab.l[i]));
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	    }
	  else{
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], (double)bab.l[j]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], (double)bab.l[i]) );
		SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	      }
	  SCIP_CALL( SCIPaddCons(scip, cons) );
	  SCIP_CALL( SCIPreleaseCons(scip, &cons) );	
	  l++;
	 }
       }


  l = cqp->ind_max_x;	
  /* y_ii >= x_i*/   
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	{
	if(i==j && ((i < cqp->nb_int) ||  (i >= cqp->nb_int && cqp->l[i]>=1)))
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
	  if(i==j && (i >= cqp->nb_int && cqp->u[i]==1 && cqp->l[i]==0))
	    {
		SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0 ));
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
		   

  
  
  /*objective function*/
  SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name,0, NULL, NULL, 0, NULL, NULL,NULL,constante,constante));
  SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, obj[0], -1.00) );
  for( i = 0; i < cqp->n; i++ )
    {
      SCIP_CALL(SCIPaddSquareCoefQuadratic(scip, cons, x[i],(cqp->new_q[ij2k(i,i,cqp->n)]- ((cqp->new_lambda_min)*2))/2));
      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], cqp->c[i] ) );
      for( j = i+1; j < cqp->n; j++ )
	{
	 SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j], cqp->new_q[ij2k(i,j,cqp->n)] )); 
      }
  }

  
  i=cqp->ind_max_x;
  for(j=0;j<cqp->ind_max_x;j++)
    for(k=j;k<cqp->ind_max_x ;k++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)])){
	     if (j==k) 
	   {
		SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i],(cqp->q[ij2k(j,k,cqp->n)] - cqp->new_q[ij2k(j,k,cqp->n)]/2 + cqp->new_lambda_min) ) );
		i++;
	   }
	 else
	   {
		SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i],2*(cqp->q[ij2k(j,k,cqp->n)] - cqp->new_q[ij2k(j,k,cqp->n)]/2 + cqp->new_lambda_min) ));
		i++;
	   }
       
       

      }
  SCIP_CALL( SCIPaddCons(scip, cons) );
  SCIP_CALL( SCIPreleaseCons(scip, &cons) );

  
  return SCIP_OKAY;
}
	       
 




int eval_sbb_miqcr_scip(struct miqcp_bab bab,double ** sol_xy, double * sol_adm, int computeLocalSol)
{
		
		
   double   sol_admissible;
   double * x_y;
   double fx;
   
   
  

   nodecount++;
   
   
   SCIP* scip;
   int i,k,l;
   
   
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, cqp->nb_col -1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
   
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   
   int feasible_for_const=1;
   
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
   
   /* Fill in the data for the problem.  */
   
   
   SCIP_CALL( setupProblem_sbb_miqcr(scip,bab));
   
   /* SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", EPS_LOCAL_SOL) ); */
   /* SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MAX_TIME_SOL_LOCAL )); */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", PRINT_LEVEL_SCIP) );

  
   // SCIP_CALL(SCIPwriteOrigProblem 	(	scip, NULL,NULL,NULL));
   
   SCIP_CALL( SCIPsolve(scip) );
   
   
   SCIP_Real * xvals;
   SCIP_CALL( SCIPallocBufferArray(scip, &xvals, cqp->nb_col -1) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),cqp->nb_col -1,x, xvals) );
   
   SCIP_Real * objval; 
   SCIP_CALL( SCIPallocBufferArray(scip, &objval, 1) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),1,obj, objval) ); 
   
   // copy the solution in our structures
   x_y = alloc_vector_d(cqp->nb_col - 1);
   for (i=0;i<cqp->nb_col - 1;i++)
     x_y[i]=xvals[i];
   
   
   sol_admissible =objval[0];
   
   
   
   
   /*******************************************************************************/
   /*******   prune branch if necessary *********************************************/
   /*******************************************************************************/
   if (sol_admissible - ub_bb  > -EPS_BB )
     {
       printf("\nnode %ld\t current sol = %6lf > UB = %lf, branch pruned ",nodecount, sol_admissible,ub_bb);
       if (nodecount==1)
	 lb_bb = sol_admissible;
       return 0;
     }
  

   if (nodecount==1)
     lb_bb = sol_admissible;
   
   int status;

   double gap;
   int is_int;
   /*Compute local sol*/
   if (computeLocalSol)
     {
       status = compute_local_sol_scip_bb(&bab.l,&bab.u);
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
   
   /* free memory arrays */
   /* for (i=0;i<cqp->nb_col - 1;i++) */
   /*   SCIP_CALL( SCIPreleaseVar(scip,&x[i]) ); */
   
   /* SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) ); */
   
   /* SCIPfreeMemoryArray(scip, &obj); */
   /* SCIPfreeMemoryArray(scip, &x); */
   
   /* SCIP_CALL( SCIPfree(&scip) ); */
   
   return 1;
}


 void sbb_miqcr_scip( struct miqcp_bab bab)
 {
		
  double   sol_admissible,fx;
  double * x_y = alloc_vector_d (cqp->nb_col - 1);
  
  int *i =  alloc_vector(1);
  int *j = alloc_vector(1);
  
  int k,m,p,ind;
  
  struct miqcp_bab new_bab_1; 
  struct miqcp_bab new_bab_2;
  father=-1;
  
  int test_time;
  int feasible_for_const;
  
  
  lb_bb = MIN_SOL_BB;
    
  /*last parameter: 1 means that it is necessary to compute a local sol at this node else 0 */
  int status = eval_sbb_miqcr_scip(bab,  &x_y,&sol_admissible,1);
  
  if(status ==1)
    {
    insert_liste(liste,sol_admissible, x_y, bab, cqp->nb_col - 1);
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
	   
	   if (test_bound_miqcp(new_bab_1,cqp->nb_col - 1))
	     {
		status = eval_sbb_miqcr_scip(new_bab_1, &x_y, &sol_admissible,1);
		if (status ==1)
		  {
		   if (ub_bb > sol_admissible)
		     {
		      absolute_gap = (sol_admissible -ub_bb);
		      relative_gap = absolute_gap/ub_bb;
		      v_abs_ref(&relative_gap);
		      if (relative_gap  > EPS_BB )
			insert_liste(liste,sol_admissible, x_y, new_bab_1, cqp->nb_col - 1);
	              }
	         }
	      }
	     
	   if (test_bound_miqcp(new_bab_2,cqp->nb_col - 1))
	     {
		status = eval_sbb_miqcr_scip(new_bab_2, &x_y, &sol_admissible,1);
		if (status ==1)
		  {
		if (ub_bb > sol_admissible)
		  {
		absolute_gap = (sol_admissible -ub_bb);
		relative_gap = absolute_gap/ub_bb;
		v_abs_ref(&relative_gap);
		if (relative_gap  > EPS_BB )
		  insert_liste(liste,sol_admissible, x_y, new_bab_2, cqp->nb_col - 1);
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
	     
		if (test_bound_miqcp(new_bab_1,cqp->nb_col - 1))
		  {
		if (x_inf_new_u)
		  status=eval_sbb_miqcr_scip(new_bab_1,  &x_y, &sol_admissible,0);
		else
		  status =eval_sbb_miqcr_scip(new_bab_1, &x_y, &sol_admissible,1);
		if (status ==1)
		  {
		if (ub_bb > sol_admissible)
		  {
		absolute_gap = (sol_admissible -ub_bb);
		relative_gap = absolute_gap/ub_bb;
		v_abs_ref(&relative_gap);
		if (relative_gap  > EPS_BB )
		  insert_liste(liste,sol_admissible, x_y, new_bab_1, cqp->nb_col - 1);
			 
	      }
	      }
	      }
	     
		if (test_bound_miqcp(new_bab_2,cqp->nb_col - 1))
		  {
		if  (x_inf_new_u)
		  status = eval_sbb_miqcr_scip(new_bab_2,  &x_y,& sol_admissible,1);
		else
		  status = eval_sbb_miqcr_scip(new_bab_2, &x_y,& sol_admissible,0);
		if (status ==1)
		  {
		if (ub_bb > sol_admissible)
		  {
		absolute_gap = (sol_admissible -ub_bb);
		relative_gap = absolute_gap/ub_bb;
		v_abs_ref(&relative_gap);
		if (relative_gap  > EPS_BB )
		  insert_liste(liste,sol_admissible, x_y, new_bab_2, cqp->nb_col - 1);
			
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


