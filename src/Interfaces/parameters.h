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

#ifndef PARAMETERS_H
#define PARAMETERS_H

/*********************************************************************************************/
/************************** Parameter for presolve with OBBT *********************************/
/*********************************************************************************************/
int OBBT_QUAD =0;     /*1 if we add quadratic constraints in the solution of Obbt problems, 0 for standard OBBT*/
int OBBT_CONT =0;     /*1 if we proceed OBBT on the constraints*/ 
int NB_PASS_OBBT_INIT= 1;  /*Max number of iteration of Obbt at the root node */
int NB_PASS_OBBT =1;  /*Max number of iteration of Obbt at each node */
double EPS_OBBT= 1e-2; /*Accuracy for updating an upper/lower bound on a variable or a left/right hand side of a constraint */
int OBBT=0;



/*********************************************************************************************/
/************************** *Define constants  ***********************************************/
/*********************************************************************************************/

double MAXU =1e6;  /* Max upper bound on a variables (used if upper bounds are not given in the data) */
double MAX_SOL_BB=1e10;  /* Initial upper bound of the problem */
double MIN_SOL_BB =-1e10; /* Initial lower bound of the problem */ 
double MAX_SOL_SDP= 1e12; /* Initial value of the sdp */ 

double EPS_BETA= 1e-6;  /*Accuracy of a Zero, positive or negative value*/



/*********************************************************************************************/
/************************** Parameter for computing upper bounds******************************/
/*********************************************************************************************/
int MAX_TIME_SOL_INIT= 5;    /*Max time for computing an upper bound at the root node */
int MAX_TIME_SOL_LOCAL= 5;   /*Max time for computing an upper bound at each node */
double EPS_LOCAL_SOL= 1e-4;  /*Accuracy for computing an upper bound */ 




/*********************************************************************************************/
/*************************** parameters for CPLEX ********************************************/
/*********************************************************************************************/

double ABS_GAP= 0.999;  /* Value of absolute gap of Cplex*/
int MEMORY =64000; /* value of CPXPARAM_Memory*/   
double OBJDIF= 0.999;  /* Value of obj_diff of Cplex*/
int TIME_LIMIT_CPLEX=3600;  /* value of CPXPARAM_TimeLimit*/   
double REL_GAP=1e-5;  /* value of CPXPARAM_MIP_Tolerances_MIPGap */    
int VARSEL=0; /*-1 or 1 or 2 or 3 or 4 default 0*/


/*********************************************************************************************/
/*************************** parameters for SCIP  ********************************************/
/*********************************************************************************************/

int PRINT_LEVEL_SCIP = 0; /*Level of print for Scip */

/*********************************************************************************************/
/*************************** parameters for IPOPT  ********************************************/
/*********************************************************************************************/

int PRINT_LEVEL_IPOPT = 0;

/*********************************************************************************************/
/*************************** parameters for CONIC BUNDLE *************************************/
/*********************************************************************************************/


double EPS_VIOL_CB= 1e-3;   /* Accuracy of the violation of a constraints*/
double EPS_TERM_CB =1e-2;   /* parameter SG_NORM of the conic bundle */
int ACT_BOUND= 1;   /* parameter Active bound of the conic bundle */
int EVAL_LIMIT= 1000;   /* parameter Eval limit of the conic bundle */
int UPDATE_LIMIT= 500;  /* parameter Update limit of the conic bundle */

double UPPER_BOUND= 1e5;
int MAX_SUBG= 1; /* Nb max subgradient added at each iteration of the SDP solver ! do not change ! */
int NB_MAX_ITER= 150; /*Max number of iteration of the SDP solver*/
int TIME_LIMIT_CB= 3600; /* Time limit for the SDP solver */
double EPS_STEP=1e-2; /* Accuracy of a step of the SDP solver */ 
double EPS_SDP=1e-4;  /* Accuracy of the SDP solver */ 
int NB_DESCENT_STEP=1; /*Number of descent steps of each iteration */
double FACTOR =1; /*Proportion of the considered constraints into the SDP solver */
int NB_MAX_CONT; /* Nb max cont considered into the sdp (default all) */

/*********************************************************************************************/
/*************************** parameters for spatial BB ***************************************/
/*********************************************************************************************/

double EPS_LM =1e-6;   /* Accuracy for considering non negative the new smallest eigenvalue */
double GAMMA =0;    /*parameter that determines the value for branching */
double EPS_INT= 1e-2;  /*Accuracy for considering that a value is integer */
double EPS_ABS_GAP= 0.99; /* Absolute gap of the spatial BB */
int NB_NODE_COMPUTE_LOCAL= 10;  /*Frequency for running local solver during the spatial BB */
double EPS_BB = 1e-5;  /*Relative gap for ending the spatial BB */
double EPS_BRANCH = 1e-2; /*Acurracy for pruning a branch in the spatial BB */
int TIME_LIMIT_BB = 3600; /*Time limit for the spatial BB */
int father;
int PRINT_LEVEL_NODE =10;
int init_param();

#endif
