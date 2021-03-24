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


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include"utilities.h"
#include"quad_prog.h"
#include "parameters.h"


int init_param()
{
 
  FILE * file_param = fopen("param.smiqp","r");
  char temp[25] = "";
  int a;
  double b;
  int ret;
  int i=0;
  
  if (file_param!=NULL)
    {
       printf("\nUse parameter file : param.smiqp\n");
      ret = fscanf(file_param, "%s", temp);
      while(i<43)
	{
	  i++;
	  if (strcmp(temp,"OBBT_QUAD")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a==0 || a==1)
	      OBBT_QUAD =a;
	  }
	  if (strcmp(temp,"OBBT_CONT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a==0 || a==1)
	      OBBT_CONT =a;
	  }
	  if (strcmp(temp,"OBBT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a==0 || a==1)
	      OBBT =a;
	  }	
	  if (strcmp(temp,"NB_PASS_OBBT_INIT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      NB_PASS_OBBT_INIT =a;
	  }
	  
	  if (strcmp(temp,"NB_PASS_OBBT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      NB_PASS_OBBT =a;
	  }
	  
	  if (strcmp(temp,"EPS_OBBT")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_OBBT =b;
	  }
	  
	  if (strcmp(temp,"MAXU")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      MAXU =b;
	  }
	  
	  if (strcmp(temp,"MAX_SOL_BB")==0){
	    ret = fscanf(file_param," %lf",&b);
	    MAX_SOL_BB=b;
	}
	  
	  if (strcmp(temp,"MIN_SOL_BB")==0){
	    ret = fscanf(file_param," %lf",&b);
	    MIN_SOL_BB=b;
	  }

	  if (strcmp(temp,"MAX_SOL_SDP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    MAX_SOL_SDP =b;
	  }
	
	  if (strcmp(temp,"EPS_BETA")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_BETA =b;
	  }
	  
	  if (strcmp(temp,"MAX_TIME_SOL_INIT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      MAX_TIME_SOL_INIT =a;
	  }
	  
	  if (strcmp(temp,"MAX_TIME_SOL_LOCAL")==0){
	    ret = fscanf(file_param," %d",&a);
	    MAX_TIME_SOL_LOCAL =a;
	  }
	  
	  if (strcmp(temp,"EPS_LOCAL_SOL")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_LOCAL_SOL =b;
	  }
	  
	  if (strcmp(temp,"ABS_GAP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      ABS_GAP =b;
	  }
	  
	  if (strcmp(temp,"MEMORY")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      MEMORY =a;
	  }
	  
	  if (strcmp(temp,"OBJDIF")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      OBJDIF =b;
	  }
	  
	  if (strcmp(temp,"TIME_LIMIT_CPLEX")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      TIME_LIMIT_CPLEX=a;
	  }
	  
	  if (strcmp(temp,"REL_GAP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      REL_GAP=b;
	  }
	  
	  
	  if (strcmp(temp,"VARSEL")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>=0 && a<=4)
	     VARSEL =a;
	  }
	  
	  if (strcmp(temp,"PRINT_LEVEL_SCIP")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>=0 && a<=5) 
	      PRINT_LEVEL_SCIP =a;
	  }
	  
	  if (strcmp(temp,"PRINT_LEVEL_IPOPT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>=0)
	      PRINT_LEVEL_IPOPT =a;
	  }
	  
	  if (strcmp(temp,"EPS_VIOL_CB")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_VIOL_CB =b;
	  }
	  
	  if (strcmp(temp,"EPS_TERM_CB")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_TERM_CB =b;
	  }
	  
	  if (strcmp(temp,"ACT_BOUND")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      ACT_BOUND =a;
	  }
	  
	  if (strcmp(temp,"EVAL_LIMIT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      EVAL_LIMIT =a;
	  }

	  if (strcmp(temp,"UPDATE_LIMIT")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      UPDATE_LIMIT =a;
	  }
	  
	  if (strcmp(temp,"UPPER_BOUND")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      UPPER_BOUND =b;
	  }
	  
	  if (strcmp(temp,"NB_MAX_ITER")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      NB_MAX_ITER=a;
	  }
	  
	  if (strcmp(temp,"TIME_LIMIT_CB")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	       TIME_LIMIT_CB =a;
	  }
	  
	  if (strcmp(temp,"EPS_STEP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_STEP =b;
	  }
	  
	  if (strcmp(temp,"EPS_SDP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	     EPS_SDP =b;
	  }
	  
	  
	  if (strcmp(temp,"FACTOR")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>=0 && b<=1)
	      FACTOR =b;
	  }
	  
	  if (strcmp(temp,"EPS_LM")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_LM =b;
	  }
	  
	  if (strcmp(temp,"GAMMA")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>=0 && b<=1)
	      GAMMA =b;
	  }
	  
	  if (strcmp(temp,"EPS_INT")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_INT =b;
	  }
	  
	  if (strcmp(temp,"EPS_ABS_GAP")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_ABS_GAP=b;
	  }
	  
	  if (strcmp(temp,"NB_NODE_COMPUTE_LOCAL")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      NB_NODE_COMPUTE_LOCAL=a;
	  }
	  
	  if (strcmp(temp,"EPS_BB")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	      EPS_BB =b;
	  }
	  
	  if (strcmp(temp,"EPS_BRANCH")==0){
	    ret = fscanf(file_param," %lf",&b);
	    if (b>0)
	    EPS_BRANCH =b;
	}
	  
	  if (strcmp(temp,"TIME_LIMIT_BB")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a >0)
	      TIME_LIMIT_BB =a;
	  }
	  
	  if (strcmp(temp,"PRINT_LEVEL_NODE")==0){
	    ret = fscanf(file_param," %d",&a);
	    if (a>0)
	      PRINT_LEVEL_NODE =a;
	  }

	  ret = fscanf(file_param, "%s", temp);
	}
      
      fclose(file_param);
      return 1;
    }
  else
    {
      printf("\nNo file for parameter found, use default parameters\n");
      return 1;
    }
 
  //fprintf(file_save,"Parameters for cplex \n TIME_LIMIT_CPLEX= %d s \n REL_GAP : %lf\n VARSEL= %d\n ABS_GAP= %lf\n OPT_TOL= %lf\n CONV_TOL= %lf\n MEMORY= %d\n BARRIER_TOL= %lf\nOBJDIF=%lf\n EPRHS= %lf\n EPRELAX= %lf\n HEURFREQ= %d\n SYMBREAK %d\n PRSOLVENODE %d\n\nParameters for CONIC BUNDLE \n PREC= %lf\n PREC_STEP= %lf\n NB_DESCENT_STEP= %d \n EPS_BETA= %lf\n EPS_LM= %lf\n EPS_TERM_CB= %lf\n ACT_BOUND= %d\n EVAL_LIMIT =%d\n UPDATE_LIMIT= %d\n FACTOR=%lf\n UPPER_BOUND=%lf\n MAX_SUBG=%d \n MAX_SUBG=%d\n NB_MAX_ITER=%d\nTIME_LIMIT_CB=%d\n Parameters for BB\n EPS_BB= %lf \n EPS_BRANCH= %lf \n TIME_LIMIT_BB= %d \n\nParameters for CSDP \nAXTOL= %lf\n AYTOL=%lf\n OBJTOL= %lf\n PINFTOL= %lf\n DINFTOL= %lf\n MAXITER= %d\n MINSF= %lf\n MAXSF= %lf\n MINSTEPP= %lf\n MINSTEPD= %lf\n USEXZGAP= %d\n TWEAKGAP= %d\n AFFINE= %d\n PERTURBOBJ= %d\n FASTMODE= %d\n\n",TIME_LIMIT_CPLEX,REL_GAP,VARSEL,ABS_GAP,OPT_TOL,CONV_TOL,MEMORY,BARRIER_TOL,OBJDIF,EPRHS,EPRELAX,HEURFREQ,SYMBREAK,PRESOLVENODE,PREC,PREC_STEP,NB_DESCENT_STEP,EPS_BETA, EPS_LM, EPS_TERM_CB,ACT_BOUND,EVAL_LIMIT,UPDATE_LIMIT,FACTOR,UPPER_BOUND,MAX_SUBG,MAX_SUBG,NB_MAX_ITER,TIME_LIMIT_CB,EPS_BB,EPS_BRANCH,TIME_LIMIT_BB,AXTOL,AYTOL,OBJTOL,PINFTOL,DINFTOL,MAXITER,MINSF,MAXSF,MINSTEPP,MINSTEPD,USEXZGAP,TWEAKGAP,AFFINE,PERTURBOBJ,FASTMODE);
       
}
