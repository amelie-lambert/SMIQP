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



#include<sys/stat.h>
#include<fcntl.h>
#include<stdio.h>
#include<stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include<string.h>
#include <time.h>
#include <math.h>

#include "solver_sdp.h"
#include "solver_sdp_mixed.h"
#include "utilities.h"
#include"quad_prog.h"

#include "cb_cinterface.h"
#include "declarations.h"
#include "solver_sdp_mixed.h"

/* parameters for CONIC BUNDLE */
extern double MAX_SOL_SDP;
extern double EPS_BETA;
extern double EPS_VIOL_CB;
extern double EPS_TERM_CB;
extern int ACT_BOUND;
extern int EVAL_LIMIT;
extern int UPDATE_LIMIT;
extern double FACTOR;
extern double UPPER_BOUND;
extern int MAX_SUBG; /* a ne pas changer sinon incohérence*/
extern int NB_MAX_ITER; /*parametre perso*/
extern int TIME_LIMIT_CB;
extern double EPS_STEP;
extern double EPS_SDP; 
extern int NB_DESCENT_STEP;


extern double sol_sdp;
extern SDP psdp;
extern long start_time;
extern int nb_var_cont; 
extern  int nb_var_int;
extern  int nb_var_01;
extern  int nb_cont_01;
extern  int nb_cont_sup_1;

void run_sdp_solver_mixed( double **x)
{
 
  /* Set parameters*/
 
    /* The problem and solution data.  */
      
  struct blockmatrix C;
  double *b;
  struct constraintmatrix *constraints;
      
  /*  Storage for the initial and final solutions. */
      
  struct blockmatrix X,Z;
  double *y;
  double pobj,dobj;
      
  /* blockptr will be used to point to blocks in constraint matrices. */
      
  struct sparseblock *blockptr;
      
  /* A return code for the call to easy_sdp().    */
  int ret;
      
  /*********************************************************************************/
  /*************************** CREATE DATA *****************************************/
  /*********************************************************************************/

  int taille;
  FILE * fsdqp;
  int cont,nb_ecart,i,j,k;
  int s;
  int r=1;
  int compt_ecart=1;
  
  int nb_cont;
  double * alpha;
  double * alphabis;
  double temp;
      
  double * q_beta;
  double * c_beta;
  double l_beta;
      
  
    
  int numnz;
  
  int file_sol;
  int stdoutsave;

  
  /********************************************************************************/
  /************************************** Dimension *******************************/
  /********************************************************************************/
 

  


  
  nb_cont= 1+ 3*(nb_var_cont +nb_var_int) + nb_var_int + nb_var_01 + nb_cont_01 + nb_cont_sup_1+psdp->m+psdp->p +psdp->mq + psdp->pq;
  nb_ecart= psdp->p+psdp->pq+ 3*(nb_var_cont +nb_var_int) + nb_var_int+ nb_cont_01 + nb_cont_sup_1;
  

  /********************************************************************************/
  /********************************** Write C matrix*******************************/
  /********************************************************************************/
  C.nblocks=2;
  C.blocks=(struct blockrec *)malloc(3*sizeof(struct blockrec));
  if (C.blocks == NULL)
    {
      printf("Couldn't allocate storage for C!\n");
      exit(1);
    };
  /*
    /   * Setup the first block.
   */
  
  C.blocks[1].blockcategory=MATRIX;
  C.blocks[1].blocksize=psdp->n +1;
  C.blocks[1].data.mat=(double *)malloc(((psdp->n+1)*(psdp->n+1))*sizeof(double));
  if (C.blocks[1].data.mat == NULL)
    {
      printf("Couldn't allocate storage for C!\n");
      exit(1);
    };
 
 
 
  dualized_objective_function(psdp, &q_beta, &c_beta, &l_beta);


  
  /*
   * Put the entries into the first block.
   */
  
  C.blocks[1].data.mat[ijtok(1,1,(psdp->n+1))]=l_beta;
  for (i=1;i<psdp->n+1;i++)
    {
      C.blocks[1].data.mat[ijtok(1,i+1,(psdp->n+1))]=c_beta[i-1]/2;
      C.blocks[1].data.mat[ijtok(i+1,1,(psdp->n+1))]=c_beta[i-1]/2;
     } 
 
  for (i=1;i<psdp->n+1;i++)
    for (j=1;j<psdp->n+1;j++)
      {
	C.blocks[1].data.mat[ijtok(i+1,j+1,(psdp->n+1))]=q_beta[ij2k(i-1,j-1,psdp->n)];
	C.blocks[1].data.mat[ijtok(j+1,i+1,(psdp->n+1))]=q_beta[ij2k(i-1,j-1,psdp->n)]	;
      }


  C.blocks[2].blockcategory=DIAG;
  C.blocks[2].blocksize=nb_ecart;
   C.blocks[2].data.vec=(double *)malloc((nb_ecart+1)*sizeof(double));
  if (C.blocks[2].data.vec == NULL)
    {
      printf("Couldn't allocate storage for C!\n");
      exit(1);
    };

  for (i=1;i<=nb_ecart;i++)
    C.blocks[2].data.vec[i]=0.0;

  


  /**********************************************************************************/
  /********************************  Constraints ************************************/
  /**********************************************************************************/
  
  /*
   * Allocate storage for the right hand side, b.
   */

  b=(double *)malloc((nb_cont+1)*sizeof(double));
  if (b==NULL)
    {
      printf("Failed to allocate storage for a!\n");
      exit(1);
    };
  
  r=1;
  
 
  /*Aq,X = bq*/
  if(!Zero_BB(psdp->mq))
    {
      for (i=0;i<psdp->mq;i++)
	{
	  b[r]=psdp->bq[i];
	  r++;
	}
    }
  
  /*DqX =eq*/
  if(!Zero_BB(psdp->pq))
    {
      for (i=0;i<psdp->pq;i++) 
	{
	  b[r]=psdp->eq[i];
	  r++;
	}
    }


  

  /* if x non binary variable */
  /*X_ii <= (u_i + l_i) x_i - u_il_i*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	b[r]=-psdp->u[i]*psdp->l[i];
	r++;
      }
  
 
  /*X_ii >= 2u_ix_i -u_i^2 */
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	b[r]=psdp->u[i]*psdp->u[i];
	r++;
      }
  
   /*X_ii >=  2l_ix_i - l_i^2 */
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	b[r]=psdp->l[i]*psdp->l[i];
	r++;
      }
  /* if x integer variable*/
  /*X_ii >= x_i */
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]!=1)
      {
	b[r]=0;
	r++;
      }
  
  /* if x continuous and l[i] >=1*/
  /*X_ii >= x_i */
  for (i=psdp->nb_int;i<psdp->n;i++)
    if (psdp->l[i]>=1)
      {
	b[r]=0;
	r++;
      }
  /* if x binary variable */
  /*X_ii = x_i */
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]==1)
      {
	b[r]=0;
	r++;
      }
  
  /* if x continuous variable between 0 an 1*/
  /*X_ii <= x_i */
  for (i=psdp->nb_int;i<psdp->n;i++)
    if ((psdp->u[i]==1) && psdp->l[i]==0)
      {
	b[r]=0;
	r++;
      }
  
  /*Ax=b*/
  if(!Zero_BB(psdp->m))
    {
      for (i=0;i<psdp->m;i++)
	{
	  b[r]=psdp->b[i];
	  r++;
	}
    }
  /*Dx<=e*/
  if(!Zero_BB(psdp->p))
    {
      for (i=0;i<psdp->p;i++)
	{
	  b[r]=psdp->e[i];
	  r++;
	}
    }
  
  b[r]=1;
  


  constraints=(struct constraintmatrix *)malloc((nb_cont+1)*sizeof(struct constraintmatrix));
  if (constraints==NULL)
    {
      printf("Failed to allocate storage for constraints!\n");
      exit(1);
    };
  
  r=1;
 
  

  /*AqX =bq*/
  /* we do not consider pure continuous products*/
  if(!Zero_BB(psdp->mq))
    for (i=0;i<psdp->mq;i++)
      {
	constraints[r].blocks=NULL;
  	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries= (psdp->n+1) * (psdp->n+2)/2;
	
  	blockptr->entries=(double *) malloc((((psdp->n+1)*(psdp->n+2))/2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((((psdp->n+1)*(psdp->n+2))/2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((((psdp->n+1)*(psdp->n+2))/2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	s=1;
  	for(j=0;j<psdp->n+1;j++)
	  for(k=j;k<psdp->n+1;k++)
	    {
	      blockptr->iindices[s]=j+1;
	      blockptr->jindices[s]=k+1;
	      blockptr->entries[s]=psdp->aq[ijk2l(i,j,k,psdp->n+1,psdp->n+1)];
	      s++;
	    }


	
	blockptr->next=constraints[r].blocks;
  	constraints[r].blocks=blockptr;

	r++;
  }


  /*x^tDqx +s = eq */
  /* we do not consider pure continuous products*/
  if(!Zero_BB(psdp->pq))
    for (i=0;i<psdp->pq;i++)
      {
  			
	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;
	
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries= ((psdp->n+1) * (psdp->n+2))/2;
	
  	blockptr->entries=(double *) malloc((((psdp->n+1) * (psdp->n+2))/2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((((psdp->n+1) * (psdp->n+2))/2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((((psdp->n+1) * (psdp->n+2))/2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
  	
	s=1;
  	for(j=0;j<psdp->n+1;j++)
	  for(k=j;k<psdp->n+1;k++)
  	  {
  	    blockptr->iindices[s]=j+1;
  	    blockptr->jindices[s]=k+1;
  	    blockptr->entries[s]=psdp->dq[ijk2l(i,j,k,psdp->n+1,psdp->n+1)];
  	    s++;
  	  }
	
        blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	r++;
      }


  /* constraints that ensure the existence of a feasible solution*/
  
  /* Xii<= uixi + lixi - uili (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	constraints[r].blocks=NULL;
  				
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	

	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=(double)-(psdp->u[i]+psdp->l[i])/2;
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;	
	
  	r++;
      }
  
  /*- Xii <=  - 2uixi + ui² (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=(double)psdp->u[i];
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=-1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;	

  	r++;
  

      }

  /* -Xii <= 2l_iu_i + l_i^2 (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=(double)psdp->l[i];
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=-1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;	

  	r++;
  

      }

  
  
  
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]!=1)
      {
  	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	constraints[r].blocks=NULL;
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=0.5 ;
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=-1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;
  	
	r++;
      }
/*- Xii <=  - xi (for integer variables of continous variable such taht l[i]>=1)*/
  for (i=psdp->nb_int;i<psdp->n;i++)
    if (psdp->l[i]>=1)
      {
  	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	constraints[r].blocks=NULL;
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=0.5 ;
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=-1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;
  	
	r++;
      }

  /*X_ii = x_i (for binary variables) */
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]==1)
      {
	constraints[r].blocks=NULL;  
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;
	
  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=-0.5;
		
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

 	r++;
      }

  /*Xii <= xi (for continuous variables beetween 0 and 1)*/
  for (i=psdp->nb_int;i<psdp->n;i++)
    if (psdp->u[i]==1 && psdp->l[i]==0)
      {
  	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	constraints[r].blocks=NULL;
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=2;

  	blockptr->entries=(double *) malloc((2+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((2+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((2+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->iindices[1]=1;
  	blockptr->jindices[1]=i+2;
  	blockptr->entries[1]=-0.5 ;
	
  	blockptr->iindices[2]=i+2;
  	blockptr->jindices[2]=i+2;
  	blockptr->entries[2]=1.0;

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;
  	
	r++;
      }

  
 /*Ax =b*/
  if(!Zero_BB(psdp->m))
    for (i=0;i<psdp->m;i++)
      {
  	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=psdp->n;
	
  	blockptr->entries=(double *) malloc((psdp->n+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((psdp->n+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((psdp->n+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	s=1;
  	for (j=0;j<psdp->n;j++)
  	  {
  	    blockptr->iindices[s]=1;
  	    blockptr->jindices[s]=j+2;
  	    blockptr->entries[s]=psdp->a[ij2k(i,j,psdp->n)]/2;
  	    s++;
  	  }

  	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	r++;

      }

  /*Dx +s = e */
  if(!Zero_BB(psdp->p))
    for (i=0;i<psdp->p;i++)
      {
  	constraints[r].blocks=NULL;
	/*
  	 * Allocate space for block 2 of constraints.
  	 */
	
	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
	blockptr->blocknum=2;
	blockptr->blocksize=nb_ecart;
	blockptr->constraintnum=r;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((1+1)*sizeof(double));
	if (blockptr->entries==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->iindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
	if (blockptr->jindices==NULL)
	  {
	    printf("Allocation of constraint block failed!\n");
	    exit(1);
	  };
	
	/*
	 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
	 */
	
	blockptr->numentries=1;
	
	/*
	 * The entry in the 1,1 position of block 3 of A1 is 1.0
	 */
	
	blockptr->iindices[1]=compt_ecart;
	blockptr->jindices[1]=compt_ecart;
	blockptr->entries[1]=1.0;
	compt_ecart++;
	
	/*
	 * Insert block 2 into the linked list of A1 blocks.  
	 */
	
	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };

	/*
  	 * Allocate space for block 1 of constraints.
  	 */
	
  	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  	if (blockptr==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
	
  	/*
  	 * Initialize block 1.
  	 */
	
  	blockptr->blocknum=1;
  	blockptr->blocksize=psdp->n+1;
  	blockptr->constraintnum=r;
  	blockptr->next=NULL;
  	blockptr->nextbyblock=NULL;
	blockptr->numentries=psdp->n;
	
  	blockptr->entries=(double *) malloc((psdp->n+1)*sizeof(double));
  	if (blockptr->entries==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  };
  	blockptr->iindices=(int *) malloc((psdp->n+1)*sizeof(int));
	
  	if (blockptr->iindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
	
  	blockptr->jindices=(int *) malloc((psdp->n+1)*sizeof(int));
  	if (blockptr->jindices==NULL)
  	  {
  	    printf("Allocation of constraint block failed!\n");
  	    exit(1);
  	  }
  	s=1;
  	for (j=0;j<psdp->n;j++)
  	  {
  	    blockptr->iindices[s]=1;
  	    blockptr->jindices[s]=j+2;
  	    blockptr->entries[s]=psdp->d[ij2k(i,j,psdp->n)]/2;
  	    s++;
  	  }

	blockptr->next=constraints[r].blocks;
	constraints[r].blocks=blockptr;

	r++;
      }
  
  
  /*X11 = 1*/
   /*
   * Allocate space for block 1 of constraints.
   */
  constraints[r].blocks=NULL;
  blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
  if (blockptr==NULL)
    {
      printf("Allocation of constraint block failed!\n");
      exit(1);
    };
  
  /*
   * Initialize block 1.
   */
  
  blockptr->blocknum=1;
  blockptr->blocksize=psdp->n+1;
  blockptr->constraintnum=r;
  blockptr->next=NULL;
  blockptr->nextbyblock=NULL;
  blockptr->numentries=1;

  blockptr->entries=(double *) malloc((1+1)*sizeof(double));
  if (blockptr->entries==NULL)
    {
      printf("Allocation of constraint block failed!\n");
      exit(1);
    };
  blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
  
  if (blockptr->iindices==NULL)
    {
      printf("Allocation of constraint block failed!\n");
      exit(1);
    }
  
  blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
  if (blockptr->jindices==NULL)
    {
      printf("Allocation of constraint block failed!\n");
      exit(1);
    }
  
  blockptr->iindices[1]=1;
  blockptr->jindices[1]=1;
  blockptr->entries[1]=1;
  
  blockptr->next=constraints[r].blocks;
  constraints[r].blocks=blockptr;
  
 


  /*************************************************************************************/
  /*************************************************************************************/
  

   
  /* Create an initial solution.  This allocates space for X, y, and Z, and sets initial values.   */
  initsoln(nb_ecart+psdp->n+1,nb_cont,C,b,constraints,&X,&y,&Z);

  ret = easy_sdp(nb_ecart+psdp->n+1,nb_cont,C,b, constraints,0.0, &X,&y,&Z, &pobj,&dobj);

  /* get back X and y values */
  double * zbeta_diag = alloc_vector_d(psdp->n);
  double * zalphaq = alloc_vector_d(psdp->mq);
  double * zalphabisq = alloc_vector_d(psdp->pq);
  
  for(i=0;i<psdp->n;i++)
    zbeta_diag[i]=0;

  s=0;
  for(i=0;i<psdp->n+1;i++)
    for(j=i;j<psdp->n+1;j++)
      {
	x[0][s]= X.blocks[1].data.mat[ijtok(i+1,j+1,(psdp->n+1))];
	s++;
      }
  

  
  s=1;
 

  /*<Aq,X> =bq > alphaq*/
  if(!Zero_BB(psdp->mq))
    for(i=0;i<psdp->mq;i++)
      {
	zalphaq[i]=y[s];
	s++;
      }
  
  /*<Dq,X> =eq -> alphabisq*/
  if(!Zero_BB(psdp->pq))
    for(i=0;i<psdp->pq;i++)
      {
	zalphabisq[i]=y[s];
	s++;
      }

  /* X_ii <= (u_i + li) x_i - uili-> beta1 (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	zbeta_diag[i]=zbeta_diag[i]+y[s];
	s++;
      }

  /* X_ii >= 2u_ix_i - u_i^2 -> beta3 (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	zbeta_diag[i] = zbeta_diag[i] - y[s];
	s++;
      }

  /* X_ii >= 2lixi - li² -> beta4 (for non binary variables)*/
  for (i=0;i<psdp->n;i++)
    if (!((psdp->u[i]==1) && (i < psdp->nb_int)))
      {
	zbeta_diag[i] = zbeta_diag[i]-y[s];
	s++;
      }

  /* X_ii >= x_i -> beta5 (for integer variables)*/
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]!=1)
      {
	zbeta_diag[i] = zbeta_diag[i]-y[s];
	s++;
      }
  
  /* X_ii >= x_i -> beta5 (for continuous variables whith l[i]>=1)*/
  for (i=psdp->nb_int;i<psdp->n;i++)
    if (psdp->l[i]>=1)
      {
	zbeta_diag[i] = zbeta_diag[i]-y[s];
	s++;
      }
  /* X_ii = x_i -> beta5  (for binary variables)*/
  for (i=0;i<psdp->nb_int;i++)
    if (psdp->u[i]==1)
      {
	zbeta_diag[i] = zbeta_diag[i] + y[s];
	s++;
      }


 /*Xii <= xi beta5(for continuous variables beetween 0 and 1)*/
  for (i=psdp->nb_int;i<psdp->n;i++)
    if (psdp->u[i]==1 && psdp->l[i]==0)
      {
	zbeta_diag[i] = zbeta_diag[i] + y[s];
	s++;
      }
	
  psdp->beta_diag=zbeta_diag;
  psdp->alphaq=zalphaq;
  psdp->alphabisq=zalphabisq;

   
}

  
int eval_fun_mixed(void* function_key, double *beta, double relprec, int max_new_subg,double *objective_value,int* new_subg,double *subgval,double *subgradient,double *x)
{
  
  /* fuction_key = SDP */
  /* beta =  variables of the problem I am solving */
  /* relprec =  relative precision of the objective values that may lead to descent steps */
  /* max_new_sub = nb of max subgradient that I had at each iteration (=1 in my case)*/
  /* objective value = value of the objective value in function of beta*/
  /* new_subg = number of new subgradients returned at each iteration (=1 in my case) */
  /* subg_values = value of the objective function for each new subgradients (=objective value in my case, because I only had 1 subgradient at my primal solution x*/
  /* subgradients = actual value of subgradients at point x (=b-B(X) in my case) */
  /* x = x variables computing during the sdp solution*/
  
  double * q_beta;
  double * c_beta;
  double l_beta;
  
  double old_objective_value=*objective_value;
  SDP psdp;
  psdp = (SDP) function_key;
  
  int dim_x= (psdp->n+2)*(psdp->n+1)/2; 
  relprec=EPS_STEP;
  max_new_subg=1;
  
  update_beta_value(psdp,beta);

  run_sdp_solver_mixed(&x);
  
  
  psdp->x=copy_vector_d(x,dim_x);
  
  dualized_objective_function(psdp, &q_beta, &c_beta, &l_beta); 
  evaluate_objective(objective_value,psdp,&q_beta,&c_beta, &l_beta);
  
 
  *new_subg=1; 

  *subgval=*objective_value;
  
  compute_subgradient(&subgradient,psdp);
 
  sol_sdp=*objective_value;
  return 0;
}




void run_conic_bundle_mixed(double * obj_val)
{
  int i;
 
  cb_problemp p;
  
  int new_length;
  
  sol_sdp= MAX_SOL_SDP;
  
  double *lower_bound=alloc_vector_d(psdp->nb_max_cont);
  double* check_violated=alloc_vector_d(psdp->length);
  int dim_x= (psdp->n+2)*(psdp->n+1)/2; 
 
  /* solution of the fisrt SDP in order to choose which constraints will be dualized*/
  double * q_beta;
  double * c_beta;
  double l;
  
  double * solution_x;
  
  int are_constraints_violated=0;

  int n_assign;
  int * assign_new_from_old=alloc_vector(psdp->nb_max_cont);
  int n_append;
  
  psdp->x=alloc_vector_d(dim_x);
  run_sdp_solver_mixed(&psdp->x);
  
  are_constraints_violated=check_constraints_violated(psdp,&check_violated,&psdp->x);
  initialize_constraints_violated(psdp,&check_violated,&new_length);
  psdp->nb_cont=new_length;
  double *x=alloc_vector_d(dim_x);			   
  double *x_center=alloc_vector_d(dim_x);			   
  double * slacks=alloc_vector_d(psdp->nb_max_cont);
  double * subgrad=alloc_vector_d(psdp->nb_max_cont);
  double *beta=alloc_vector_d(psdp->nb_max_cont);

  int test_time;
  
  /* Interface with the conic bundle method*/  
  p = cb_construct_problem(1);
  cb_init_problem(p,psdp->nb_cont,lower_bound,0);
  cb_add_function(p,(void *)psdp,(cb_functionp)eval_fun_mixed,(cb_subgextp)gen_subg,dim_x);


  free_vector_d(lower_bound);
  
  /* parameters  */
  /*precision de terminaison*/
  cb_set_term_relprec(p, EPS_SDP);
  cb_set_print_level(p,1);
  /*Nb max d evaluation de la fonction */
  cb_set_eval_limit(p,EVAL_LIMIT);
  /*Fixe certaines variables si beta tres grand*/
  cb_set_active_bounds_fixing(p, ACT_BOUND);
  /*Nb max de mise a jour si null step */
  cb_set_inner_update_limit(p,UPDATE_LIMIT); 		

  int current_dim;
  int k=1;
  /*stop = 0;*/

  add_constraints_violated(psdp,&check_violated,&n_append);

  if(n_append>0)
      {
	lower_bound=alloc_vector_d(n_append);
    	cb_append_variables(p,n_append,lower_bound,NULL);
    	free_vector_d(lower_bound);
      }

  cb_do_descent_step(p);

  int ret =0; 
  do{
    for(i=0;i<NB_DESCENT_STEP;i++)
      ret = cb_do_descent_step(p);

    if (ret==0)
      {

	cb_get_center(p,beta);
	cb_get_approximate_primal(p,(void *)psdp,x);
	cb_get_center_primal(p,(void *)psdp,x_center);
	current_dim=cb_get_dim(p);
	cb_get_approximate_slacks(p,slacks);
	
	/* Keep constraints still violated (i.e. when slack[i]==0) */
	
	purge_constraints_not_violated(psdp, &slacks, &n_assign, &assign_new_from_old);
	cb_reassign_variables(p,n_assign,assign_new_from_old);
	
	/*Add new constraints that are violated */
	are_constraints_violated=check_constraints_violated(psdp,&check_violated, &psdp->x);
	
	
	add_constraints_violated(psdp,&check_violated,&n_append);
	
	if(n_append>0)
	  {
	    lower_bound=alloc_vector_d(n_append);
	    cb_append_variables(p,n_append,lower_bound,NULL);
	    free_vector_d(lower_bound);
	  }
	
	*obj_val = cb_get_objval(p);
	k++;
      }
	test_time = time(NULL);
      
  } while ( !(cb_termination_code(p)) && (cb_get_sgnorm(p) > EPS_TERM_CB) && (k < NB_MAX_ITER) && (test_time - start_time < TIME_LIMIT_CB));
  
  
  cb_get_center(p,beta);
  update_beta_value(psdp,beta);
  run_sdp_solver_mixed(&x);
  psdp->x=copy_vector_d(x,dim_x);
    
   
  cb_print_termination_code(p);
  cb_destruct_problem(&p);		  

}


void  compute_alpha_beta_sbb_miqcr( ) 
{
 
  double obj_val;

  int i,j;

  double  *beta_bis;
  double  *zbeta;
  double * zalphaq;
  double * zalphabisq;
  double zalpha=0;
  double zalphabis=0;
   
  printf("\nStart running SDP solver\n\n");   
  zbeta = alloc_matrix_d(psdp->n,psdp->n);
  zalphaq= alloc_vector_d(psdp->mq);
  zalphabisq= alloc_vector_d(psdp->pq);
  
  
  run_conic_bundle_mixed(&obj_val);
  
  printf("\nObjective function value sdp : %lf\n\n",-obj_val);
            
  for(i=0;i<psdp->n;i++)
    for(j=i;j<psdp->n;j++)
      if (i==j) 
	zbeta[ij2k(i,i,psdp->n)]= psdp->beta_diag[i];
      else
	zbeta[ij2k(i,j,psdp->n)]= psdp->beta_1[ij2k(i,j,psdp->n)] + psdp->beta_2[ij2k(i,j,psdp->n)] - psdp->beta_3[ij2k(i,j,psdp->n)] - psdp->beta_4[ij2k(i,j,psdp->n)];
  
  for(i=0;i<psdp->n;i++)
    for(j=i+1;j<psdp->n;j++)
      zbeta[ij2k(j,i,psdp->n)]=zbeta[ij2k(i,j,psdp->n)];
      
  for(i=0;i<psdp->n ;i++)
    for(j=0;j<psdp->n ;j++)
      if (((zbeta[ij2k(i,j,psdp->n )] < EPS_BETA) && (zbeta[ij2k(i,j,psdp->n )] > 0)) || ((zbeta[ij2k(i,j,psdp->n )] > -EPS_BETA) && (zbeta[ij2k(i,j,psdp->n )] < 0)))
	zbeta[ij2k(i,j,psdp->n )]=0;
  
      zalphaq=psdp->alphaq;
      zalphabisq=psdp->alphabisq;
      
      psdp->beta = zbeta;
  
      psdp->alphaq=zalphaq;
      psdp->alphabisq=zalphabisq;

}

