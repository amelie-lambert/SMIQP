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
#include"boxqp_scip.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

extern C_MIQCP cqp;


extern long nodecount;

extern double * best_x;
extern double ub_bb;

extern int TIME_LIMIT_BB;
extern int PRINT_LEVEL_SCIP;

extern double EPS_BB;
extern double EPS_BETA;
extern double EPS_ABS_GAP;

extern int IS_CPLEX;

SCIP_VAR ** x;
SCIP_VAR ** obj;
extern double * best_x;
extern double ub_bb;

static SCIP_RETCODE setupProblem_boxqp(SCIP* scip)
{
  char name[SCIP_MAXSTRLEN];
  SCIP_CONS* cons;
  SCIP_Real constante = -cqp->cons;
    
    
 
  int i, j, k, l;
 
  /* create empty problem */
  SCIP_CALL( SCIPcreateProbBasic(scip, "boxqp") );
 
  /* change to maximization if optimizing number of circles instead of rectangle area */
  SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   
  
  /*Allocation memoire du vecteur de variable*/
  
  SCIP_CALL( SCIPallocMemoryArray(scip, &x, cqp->nb_col -1) );
    
  for( i = 0; i <cqp->nb_int; i++ )
    {
      
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, cqp->l[i], cqp->u[i], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }
    
  for( i; i <cqp->nb_col -1; i++ )
    {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, cqp->l[i], cqp->u[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }
    
  /*variable obj for the objective function*/
  SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
  SCIP_CALL( SCIPcreateVarBasic(scip, &obj[0], "objective",-SCIPinfinity(scip) ,SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
  SCIP_CALL( SCIPaddVar(scip, obj[0]) );
  
  
 
  
  l = cqp->ind_max_x;
  /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      {
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	  {
	    if (Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	    {
	      if(i==j)
		  l++;
		
	      else
		{
		  SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0));
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], -1 ));
		
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
		  SCIP_CALL( SCIPaddCons(scip, cons) );
		  SCIP_CALL( SCIPreleaseCons(scip, &cons) ); 
		  l++;
		}
	  }
	else
	  l++;
	  }
      }
      
		  
  l = cqp->ind_max_x;	    
  /* y_ij <= l_jx_i + u_ix_j - l_ju_i*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	{
	  if (Negative_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	    {
	      if (i==j)
		l++;
	      else
		{
		  SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0));
		  
		
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], -1));
		  SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], 1 ));
		  SCIP_CALL( SCIPaddCons(scip, cons) );
		  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
		  l++;
		}
	    }
	  else
	    l++;
	}
		

  l = cqp->ind_max_x;	  
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      {
       if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	 {
	   if (Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)])) 
	     {
	       if(i==j)
		 l++;
	       else
		 {
		   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 1));
		   SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 1 ));
		   SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[j], 1));
		   SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
		   SCIP_CALL( SCIPaddCons(scip, cons) );
		   SCIP_CALL( SCIPreleaseCons(scip, &cons) ); 
		   l++;
		 }
	     }
	   else
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
	   if (Positive_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	     {
	       if(i==j)
		 l++;
	       else
		 {
		   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),0));
		  
		   SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
		   SCIP_CALL( SCIPaddCons(scip, cons) );
		   SCIP_CALL( SCIPreleaseCons(scip, &cons) );	
		   l++;
		 }
	     }
	   else
	     l++;
	 }
      }


  l = cqp->ind_max_x;	
  /* y_ii >= x_i*/   
  for(i=0;i<cqp->ind_max_x;i++)
    for(j=i;j<cqp->ind_max_x;j++)
      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,cqp->n)]))
	{
	  if(i==j)
	    {
	      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0,0 ));
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[i], 1));
	      SCIP_CALL(  SCIPaddCoefLinear(scip, cons, x[l], -1 ));
	      l++;
	      SCIP_CALL( SCIPaddCons(scip, cons) );
	      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	    }
	  else
	    l++;
	  
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
      if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,cqp->n)]))
	{
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
	       
 




int boxqp_scip()
{
		
		

   SCIP* scip;
   int i,k,l;
   
   
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, cqp->nb_col -1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
   
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   
  
   
   /* Fill in the data for the problem.  */
   
   
   SCIP_CALL( setupProblem_boxqp(scip));
   PRINT_LEVEL_SCIP = 3;
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", EPS_BB) ); 
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", TIME_LIMIT_BB )); 
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel",PRINT_LEVEL_SCIP) );

   
   SCIP_CALL( SCIPsolve(scip) );
   
   
   SCIP_Real * xvals;
   SCIP_CALL( SCIPallocBufferArray(scip, &xvals, cqp->nb_col -1) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),cqp->nb_col -1,x, xvals) );
   
   SCIP_Real * objval; 
   SCIP_CALL( SCIPallocBufferArray(scip, &objval, 1) );
   SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),1,obj, objval) ); 

   nodecount =  SCIPgetNNodes(scip);
 
   // copy the solution in our structures
   for (i=0;i<cqp->nb_col -1;i++)
     best_x[i]=xvals[i];
   ub_bb =objval[0] + cqp->cons;
   
   
   /* free memory arrays */
   for (i=0;i<cqp->nb_col -1;i++)
     SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
   
   SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
   
   SCIPfreeMemoryArray(scip, &obj);
   SCIPfreeMemoryArray(scip, &x);
   
   SCIP_CALL( SCIPfree(&scip) );
   
   return 1;
}
