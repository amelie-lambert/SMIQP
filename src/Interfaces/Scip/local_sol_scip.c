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
#include "utilities.h"
#include"quad_prog.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "local_sol_scip.h"


extern int MAX_TIME_SOL_INIT;
extern int MAX_TIME_SOL_LOCAL;
extern double EPS_LOCAL_SOL;
extern int PRINT_LEVEL_SCIP;

extern MIQCP qp;
extern C_MIQCP cqp;
extern int is_not_0_1;
extern int type;


SCIP_VAR ** x;
SCIP_VAR ** obj;

/* evaluation of an initial upper bound*/
static SCIP_RETCODE setupProblem_init(SCIP* scip)
 {
    char name[SCIP_MAXSTRLEN];
    SCIP_CONS* cons;
    SCIP_Real constante = -qp->cons;

  
 
    int i, j, k;
 
    /* create empty problem */
    SCIP_CALL( SCIPcreateProbBasic(scip, "local sol init") );
 
    /* change to maximization if optimizing number of circles instead of rectangle area */
    SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   
 
    /*Allocation memoire du vecteur de variable*/
    
    SCIP_CALL( SCIPallocMemoryArray(scip, &x, qp->n) );
    
    for( i = 0; i <qp->nb_int; i++ )
    {
    
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, qp->l[i], qp->u[i], 0.0, SCIP_VARTYPE_INTEGER) );
	SCIP_CALL( SCIPaddVar(scip, x[i]) );
      }
    
    for( i; i <qp->n; i++ )
    {
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, qp->l[i], qp->u[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
       SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }

    /*variable obj for the objective function*/
    SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
    SCIP_CALL( SCIPcreateVarBasic(scip, &obj[0], "objective",-SCIPinfinity(scip) ,SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
    SCIP_CALL( SCIPaddVar(scip, obj[0]) );
    
 
    
    /* linear constraint: sum_i Ax_i = b  */
    for (j=0;j<qp->m;j++){
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name,0, NULL, NULL, qp->b[j],qp->b[j]) );
      for (i=0;i<qp->n;i++)
	if (!Zero_BB(qp->a[ij2k(j,i,qp->n)])){
	  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qp->a[ij2k(j,i,qp->n)] ) );
	}
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
    
       /* linear constraint: sum_i Dx_i <= e  */
       for (j=0;j<qp->p;j++){
	 SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),qp->e[j]) );
	 for (i=0;i<qp->n;i++)
	   if (!Zero_BB(qp->d[ij2k(j,i,qp->n)])){
	     SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qp->d[ij2k(j,i,qp->n)] ) );
	     
	   }
	 SCIP_CALL( SCIPaddCons(scip, cons) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	 
       }


       /*Quadratic constraints*/
     
       /* x^TAqx = bq  */
       for (k=0;k<qp->mq;k++)
	 {
	   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name, 0, NULL, NULL, 0, NULL, NULL,NULL, qp->bq[k], qp->bq[k]) );
	   
	   for( i = 0; i < qp->n;i++ )
	     {
	       if (!Zero_BB(qp->aq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->aq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)] ) );
	       if (!Zero_BB(qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], 2*qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)] ) );
	       for( j = i+1; j < qp->n; j++ )
		 {
		   if (!Zero_BB(qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]))
		     SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] ) ); 
		 }
	     } 
	   SCIP_CALL( SCIPaddCons(scip, cons) );
	   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
	 }

       /* x^TDqx <= eq  */
       for (k=0;k<qp->pq;k++)
	 {
	   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name, 0, NULL, NULL, 0, NULL, NULL,NULL,-SCIPinfinity(scip) , qp->eq[k]) );
	   
	   for( i = 0; i < qp->n; i++ )
	     {
	       if (!Zero_BB(qp->dq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->dq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)] ) );
	       if (!Zero_BB(qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], 2*qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)] ) );
	       for( j = i+1; j < qp->n; j++ )
		 {
		   if (!Zero_BB(qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]))
		     SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] ) ); 
		 }
	     }
	   
	   SCIP_CALL( SCIPaddCons(scip, cons) );
	   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
	 }


       /*objective function*/
       SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name,0, NULL, NULL, 0, NULL, NULL,NULL,constante,constante));
       SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, obj[0], -1.00) );
       for( i = 0; i < qp->n; i++ )
	 {
	   if (!Zero_BB(qp->q[ij2k(i,i,qp->n)]))
	     SCIP_CALL(SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->q[ij2k(i,i,qp->n)]) );
	   if (!Zero_BB(qp->c[i]))
	     SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], qp->c[i] ) );
	   for( j = i+1; j < qp->n; j++ )
	     {
	       if (!Zero_BB(qp->q[ij2k(i,j,qp->n)]))
		 SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->q[ij2k(i,j,qp->n)]) ); 
	     }
	 }
	       
       SCIP_CALL( SCIPaddCons(scip, cons) );


      
       
       SCIP_CALL( SCIPreleaseCons(scip, &cons) );
       
       
    return SCIP_OKAY;
 }
 
 
int compute_local_sol_scip_init ()
 {
    SCIP* scip;
    int i;
 
   
    
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    
 
    SCIP_CALL( setupProblem_init(scip) );
 
    SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", EPS_LOCAL_SOL) );
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MAX_TIME_SOL_INIT ));
    SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel",PRINT_LEVEL_SCIP ) );

    SCIPwriteOrigProblem 	( scip, NULL,NULL,1); 		
   SCIP_CALL( SCIPsolve(scip) );
    
    
    qp->local_sol = alloc_vector_d(qp->n);

    SCIP_Real * xvals;
    SCIP_CALL( SCIPallocBufferArray(scip, &xvals, qp->n) );
    SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),qp->n,x, xvals) );
    
    if (xvals ==NULL)
      {
	
	/* free memory arrays */
	for (i=0;i<qp->n;i++)
	  SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
	
	SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
	
	SCIPfreeMemoryArray(scip, &obj);
	SCIPfreeMemoryArray(scip, &x);
	
	SCIP_CALL( SCIPfree(&scip) );
	return 0;
      }

     SCIP_Real * objval; 
    SCIP_CALL( SCIPallocBufferArray(scip, &objval, 1) );
    SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),1,obj, objval) ); 
    
    
    for (i=0;i<qp->n;i++)
      qp->local_sol[i]=xvals[i];
    
    qp->sol_adm=objval[0];


         
    /* free memory arrays */
    for (i=0;i<qp->n;i++)
      SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
    
  SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
  
  SCIPfreeMemoryArray(scip, &obj);
  SCIPfreeMemoryArray(scip, &x);
  
  SCIP_CALL( SCIPfree(&scip) );
  
  return 1;
  
 }
 
/* evaluation of an upper bound during the spatial bb*/
static SCIP_RETCODE setupProblem_bb(SCIP* scip,double ** l, double ** u)
 {
    char name[SCIP_MAXSTRLEN];
    SCIP_CONS* cons;
    SCIP_Real constante = -qp->cons;
    
    
 
    int i, j, k;
 
    /* create empty problem */
    SCIP_CALL( SCIPcreateProbBasic(scip, "miqcr_bb") );
 
    /* change to maximization if optimizing number of circles instead of rectangle area */
    SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   
 
    /*Allocation memoire du vecteur de variable*/
    
    SCIP_CALL( SCIPallocMemoryArray(scip, &x, qp->n) );
    
    for( i = 0; i <qp->nb_int; i++ )
    {
    
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, l[0][i], u[0][i], 0.0, SCIP_VARTYPE_INTEGER) );
	SCIP_CALL( SCIPaddVar(scip, x[i]) );
      }
    
    for( i; i <qp->n; i++ )
    {
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, l[0][i], u[0][i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
       SCIP_CALL( SCIPaddVar(scip, x[i]) );
    }

    /*variable obj for the objective function*/
    SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
    SCIP_CALL( SCIPcreateVarBasic(scip, &obj[0], "objective",-SCIPinfinity(scip) ,SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
    SCIP_CALL( SCIPaddVar(scip, obj[0]) );

  
 
    
     
       /* linear constraint: sum_i Ax_i = b  */
       for (j=0;j<qp->m;j++){
	 SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name,0, NULL, NULL, qp->b[j],qp->b[j]) );
	 for (i=0;i<qp->n;i++)
	   if (!Zero_BB(qp->a[ij2k(j,i,qp->n)])){
	       SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qp->a[ij2k(j,i,qp->n)] ) );
	   }
	 SCIP_CALL( SCIPaddCons(scip, cons) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cons) );
       }

       /* linear constraint: sum_i Dx_i <= e  */
       for (j=0;j<qp->p;j++){
	 SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip),qp->e[j]) );
	 for (i=0;i<qp->n;i++)
	   if (!Zero_BB(qp->d[ij2k(j,i,qp->n)])){
	     SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],qp->d[ij2k(j,i,qp->n)] ) );

	   }
	 SCIP_CALL( SCIPaddCons(scip, cons) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	 
       }


       /*Quadratic constraints*/
        /* x^TAqx = bq  */
       for (k=0;k<qp->mq;k++)
	 {
	   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name, 0, NULL, NULL, 0, NULL, NULL,NULL, qp->bq[k], qp->bq[k]) );
	   
	   for( i = 0; i < qp->n;i++ )
	     {
	       if (!Zero_BB(qp->aq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->aq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)] ) );
	       if (!Zero_BB(qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], 2*qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)] ) );
	       for( j = i+1; j < qp->n; j++ )
		 {
		   if (!Zero_BB(qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]))
		     SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] ) ); 
		 }
	     } 
	   SCIP_CALL( SCIPaddCons(scip, cons) );
	   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
	 }

       /* x^TDqx <= eq  */
       for (k=0;k<qp->pq;k++)
	 {
	   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name, 0, NULL, NULL, 0, NULL, NULL,NULL,-SCIPinfinity(scip) , qp->eq[k]) );
	   
	   for( i = 0; i < qp->n; i++ )
	     {
	       if (!Zero_BB(qp->dq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->dq[ijk2l(k,i+1,i+1,qp->n+1,qp->n+1)] ) );
	       if (!Zero_BB(qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]))
		 SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], 2*qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)] ) );
	       for( j = i+1; j < qp->n; j++ )
		 {
		   if (!Zero_BB(qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]))
		     SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] ) ); 
		 }
	     }
	   
	   SCIP_CALL( SCIPaddCons(scip, cons) );
	   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	   
	 }


       /*objective function*/
       SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name,0, NULL, NULL, 0, NULL, NULL,NULL,constante,constante));
       SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, obj[0], -1.00) );
       for( i = 0; i < qp->n; i++ )
	 {
	   if (!Zero_BB(qp->q[ij2k(i,i,qp->n)]))
	     SCIP_CALL(SCIPaddSquareCoefQuadratic(scip, cons, x[i],qp->q[ij2k(i,i,qp->n)]) );
	   if (!Zero_BB(qp->c[i]))
	     SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, x[i], qp->c[i] ) );
	   for( j = i+1; j < qp->n; j++ )
	     {
	       if (!Zero_BB(qp->q[ij2k(i,j,qp->n)]))
		 SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j],2*qp->q[ij2k(i,j,qp->n)]) ); 
	     }
	 }
	       
       SCIP_CALL( SCIPaddCons(scip, cons) );

       SCIP_CALL( SCIPreleaseCons(scip, &cons) );
	 
  

       
       
    return SCIP_OKAY;
 }
 
 
int compute_local_sol_scip_bb (double ** l, double ** u)
 {
    SCIP* scip;
    int i;
 
  
    
    SCIP_CALL( SCIPallocMemoryArray(scip, &x, qp->n) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &obj, 1) );
    
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    
 
   
    SCIP_CALL( setupProblem_bb(scip,l, u) );
  
    SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", EPS_LOCAL_SOL) );
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", MAX_TIME_SOL_LOCAL ));
    SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", PRINT_LEVEL_SCIP) );

	       
    SCIP_CALL( SCIPsolve(scip) );
    
    
    qp->local_sol = alloc_vector_d(qp->n);
    SCIP_Real * xvals;
    SCIP_CALL( SCIPallocBufferArray(scip, &xvals, qp->n) );
    SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),qp->n,x, xvals) );
    
    if (xvals ==NULL)
      {
	
	/* free memory arrays */
	for (i=0;i<qp->n;i++)
	  SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
	
	SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
	
	SCIPfreeMemoryArray(scip, &obj);
	SCIPfreeMemoryArray(scip, &x);
	
	SCIP_CALL( SCIPfree(&scip) );
	return 0;
      }

      
   
    
    SCIP_Real * objval; 
    SCIP_CALL( SCIPallocBufferArray(scip, &objval, 1) );
    SCIP_CALL( SCIPgetSolVals(scip, SCIPgetBestSol(scip),1,obj, objval) ); 
    
    
    for (i=0;i<qp->n;i++)
      qp->local_sol[i]=xvals[i];
    
    qp->sol_adm=objval[0];

         
    /* free memory arrays */
    for (i=0;i<qp->n;i++)
      SCIP_CALL( SCIPreleaseVar(scip,&x[i]) );
    
  SCIP_CALL( SCIPreleaseVar(scip,&obj[0]) );
  
  SCIPfreeMemoryArray(scip, &obj);
  SCIPfreeMemoryArray(scip, &x);
  
  SCIP_CALL( SCIPfree(&scip) );
  
  return 1;
  
 }
 
