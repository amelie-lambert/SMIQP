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

#include"quad_prog.h"
#include "utilities.h"

#include<stdlib.h>
#include<math.h>
#include <stdio.h>
#include <string.h>

extern C_MIQCP cqp;
extern SDP psdp;
extern MIQCP qp;
extern Q_MIQCP qqp;
extern double FACTOR;
extern double MAX_SOL_BB;
extern double EPS_BETA;

extern int is_not_0_1;
extern int nb_var_cont; 
extern  int nb_var_int;
extern  int nb_var_01;
extern  int nb_cont_01;
extern  int nb_cont_sup_1;

/******************************************************************************/
/******************* Structure for classical MIQCP ***************************/
/******************************************************************************/

MIQCP new_miqcp()
{
  return (MIQCP) malloc(sizeof(struct miqcp));
}

Q_MIQCP new_q_miqcp()
{
  return (Q_MIQCP) malloc(sizeof(struct q_miqcp));
}

void update_qp()
{
  if (cqp->p <qp->p)
    {
      qp->p = cqp->p;
      free(qp->d);
      free(qp->e);
      qp->d=copy_matrix_d(cqp->d,cqp->p,cqp->n);
      qp->e=copy_vector_d(cqp->e,cqp->p);
    }

  if (cqp->pq <qp->pq)
    {
      qp->pq = cqp->pq;
      free(qp->dq);
      free(qp->eq);
      qp->dq=copy_matrix_3d(cqp->dq,cqp->pq,cqp->n+1,cqp->n+1);
      qp->eq=copy_vector_d(cqp->eq,cqp->pq);
      
    }

  free(qp->u);
  free(qp->l);
  qp->u=copy_vector_d(cqp->u,cqp->n + cqp->n *(cqp->n+1)/2);
  qp->l=copy_vector_d(cqp->l,cqp->n + cqp->n *(cqp->n+1)/2);

}




/*******************************************************************************/
/******************* Structure for general MIQCP ***************************/
/******************************************************************************/


void create_q_miqcp( )
{
  int i,j,k,l,lim;
  int test =0;
  
  qqp = new_q_miqcp();
  qqp->n=qp->n;
  qqp->nb_int=qp->nb_int;
  qqp->m=qp->m;
  qqp->p=qp->p;
  qqp->mq=qp->mq;
  qqp->pq=qp->pq;
  qqp->u = copy_vector_d(qp->u,qp->n + (qp->n+1)*qp->n/2);
  qqp->l = copy_vector_d(qp->l,qp->n + (qp->n+1)*qp->n/2);
  
  qqp->cons = qp->cons;
  qqp->q = copy_matrix_d(qp->q, qp->n,qp->n);
  qqp->c = copy_vector_d(qp->c,qp->n);
  
  qp->N = vector_length_base_2(qp->u, qp->nb_int);
  qp->lgc = create_u_base2_cumulated(qp->u,qp->n);
  qp->lg = create_u_base2(qp->u,qp->n);

  if(!Zero_BB(qp->m))
    {
      qqp->a = copy_matrix_d(qp->a, qp->m,qp->n);
      qqp->b = copy_vector_d(qp->b,qp->m);
    }  
  else
    {
      qqp->a=alloc_matrix_d(1,1);
      qqp->b = copy_vector_d(qp->b,qp->m);
    }
  
  if(!Zero_BB(qp->p))
    {
      qqp->d= copy_matrix_d(qp->d, qp->p,qp->n);
      qqp->e = copy_vector_d(qp->e,qp->p);
    }
  else
    {
      qqp->d=alloc_matrix_d(1,1);
      qqp->e = copy_vector_d(qp->e,qp->p);
    }
  
  if(!Zero_BB(qp->mq))
    {
      qqp->aq = copy_matrix_3d(qp->aq, qp->mq,qp->n+1,qp->n+1);
      qqp->bq = copy_vector_d(qp->bq,qp->mq);
    }  
  else
    {
      qqp->aq=alloc_matrix_d(1,1);
      qqp->bq = copy_vector_d(qp->bq,qp->mq);
    }
  
  if(!Zero_BB(qp->pq))
    {
      qqp->dq= copy_matrix_3d(qp->dq, qp->pq,qp->n+1,qp->n+1);
      qqp->eq = copy_vector_d(qp->eq,qp->pq);
    }
  else
    {
      qqp->dq=alloc_matrix_d(1,1);
      qqp->eq = copy_vector_d(qp->eq,qp->pq);
    }
  
  if(!Zero_BB(qp->mq))
    qqp->alphaq = copy_vector_d(psdp->alphaq,qp->mq);
  else
    qqp->alphaq = alloc_vector_d(1);

  if(!Zero_BB(qp->pq))
    qqp->alphabisq =copy_vector_d(psdp->alphabisq,qp->pq);
  else
    qqp->alphabisq = alloc_vector_d(1);
  
  qqp->beta = copy_matrix_d(psdp->beta,qp->n,qp->n);
   
  int cont = 0;
  for(i=qp->n;i>qp->n - qp->nb_int;i--)
    cont=cont+i;
  
  int cont_bis = 0;
  for(i=qp->n-1;i>qp->n- 1 - qp->nb_int;i--)
    cont_bis=cont_bis+i;
  
  /*Computation of the initial number of constraints */
  /* Ax=b */
  qqp->ind_eg = qp->m;
  /* Dx = e */
  qqp->ind_ineg = qqp->ind_eg + qp->p;
  /*x^TAqx = bq*/
  qqp->ind_eg_q = qqp->ind_ineg + qp->mq;
   /*x^TDqx = eq*/
  qqp->ind_ineg_q = qqp->ind_eg_q + qp->pq;
  /* x = sum_k 2^k t_ik*/
  qqp->ind_dec_x=qqp->ind_ineg_q+qp->nb_int;
  /* y_ij = sum_k 2^k z_ijk*/
  qqp->ind_dec_y=qqp->ind_dec_x + qp->nb_int*(qp->n);
  /*z_ijk <= u_jt_ik*/
  qqp->ind_z_1=qqp->ind_dec_y +(qp->n)*qp->N;
  /*z_ijk <= x_j*/
  qqp->ind_z_2=qqp->ind_z_1+(qp->n)*qp->N;
  /*z_ijk >= x_j - u_j(1 - t_ik)*/
  qqp->ind_z_3=qqp->ind_z_2+(qp->n)*qp->N;
  /*z_ijk >= 0*/
  qqp->ind_z_4=qqp->ind_z_3+(qp->n)*qp->N;
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j */
  qqp->ind_y_1=qqp->ind_z_4 +cont;
  /* y_ii >= x_i */
  qqp->ind_y_2= qqp->ind_y_1 + qp->nb_int;
  /* y_ij = y_ji */
  qqp->ind_y_3 =qqp->ind_y_2 + cont_bis;
  /* y_ij <= u_ix_j  for x_i integer and x_j real*/
  qqp->nrow =qqp->ind_y_3 + qp->nb_int * (qp->n - qp->nb_int);
 
  /* initialisation of the new number of constraints*/
  /* Ax=b */
  qqp->nb_eg =qp->m;
  /* Dx = e */
  qqp->nb_ineg = qp->p;
  /*x^TAqx = bq*/
  qqp->nb_eg_q = qp->mq;
  /*x^TDqx = eq*/
  qqp->nb_ineg_q = qp->pq;
   /* x = sum_k 2^k t_ik*/
  qqp->nb_dec_x=qp->nb_int;
  /* y_ij = sum_k 2^k z_ijk*/
  qqp->nb_dec_y=0;
  /*z_ijk <= u_jt_ik*/
  qqp->nb_z_1=0;
  /*z_ijk <= x_j*/
  qqp->nb_z_2=0;
  /*z_ijk >= x_j - u_j(1 - t_ik)*/
  qqp->nb_z_3=0;
  /*z_ijk >= 0*/
  qqp->nb_z_4=0;
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j */
  qqp->nb_y_1=0;
  /* y_ii >= x_i */
  qqp->nb_y_2=0;
  /* y_ij = y_ji */
  qqp->nb_y_3 =0;
  /* y_ij <= u_ix_j  for x_i integer and x_j real*/
  qqp->nb_y_4=0;

  /* real number of variables y and z*/
  qqp->nb_var_y = alloc_matrix(qp->n,qp->n);
 
  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      {
	/*initialisation*/
	qqp->nb_var_y[ij2k(i,j,qp->n)] =0;
	if (!Zero_BB(qqp->beta[ij2k(i,j,qp->n)]))
	  qqp->nb_var_y[ij2k(i,j,qp->n)] =1;
      }
  
  for(l=0;l<qp->mq;l++)
    for(i=0;i<qp->n;i++)
      for(j=0;j<qp->n;j++)
	if (!Zero_BB(qp->aq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)]) && ((i < qqp->nb_int) || (j < qqp->nb_int)))
	  qqp->nb_var_y[ij2k(i,j,qp->n)] =1;
  
  for(l=0;l<qp->pq;l++)
    for(i=0;i<qp->n;i++)
      for(j=0;j<qp->n;j++)
	if (!Zero_BB(qp->dq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)]) && ((i < qqp->nb_int) || (j < qqp->nb_int)))
	  qqp->nb_var_y[ij2k(i,j,qp->n)] =1;
  
  qqp->nb_y = nb_non_zero_matrix_i(qqp->nb_var_y, qp->n, qp->n);
  
  qqp->nb_z=0;
      for(i=0;i<qp->nb_int;i++)
	{
	  lim=log_base_2(qp->u[i])+1;
	  for(j=0;j<qp->n;j++)
	    if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	      qqp->nb_z=qqp->nb_z+lim;
	}
      
 
  /*computation of the new number of constraints */
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for (j=i;j<qp->n;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->nb_dec_y++;
	    qqp->nb_y_1++;
	    qqp->nb_z_1+=lim;
	    qqp->nb_z_2+=lim;
	    qqp->nb_z_3+=lim;
	    qqp->nb_z_4+=lim;
	    if (j>=qp->nb_int)
	      qqp->nb_y_4++;
	    if (i==j)
	      qqp->nb_y_2++;
	    else
	      qqp->nb_y_3++;
	  }
    }

  /* z variables from the triangle inf of the matrix*/
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for (j=0;j<i;j++)
	if ((!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)])) && (i!=j))
	  {
	    qqp->nb_dec_y++;
	    qqp->nb_z_1+=lim;
	    qqp->nb_z_2+=lim;
	    qqp->nb_z_3+=lim;
	    qqp->nb_z_4+=lim;
	  }
    }
 
  qqp->nb_row =qqp->nb_eg + qqp->nb_ineg + qqp->nb_eg_q + qqp->nb_ineg_q + qqp->nb_dec_x + qqp->nb_dec_y + qqp->nb_z_1+ qqp->nb_z_2+ qqp->nb_z_3+ qqp->nb_z_4+  qqp->nb_y_1+  qqp->nb_y_2 + qqp->nb_y_3+ qqp->nb_y_4;
  
 qqp->cont = alloc_vector(qqp->nrow);
 for(i=0;i<qqp->nrow;i++)
   qqp->cont[i]=-1;
 

  /* number of initial variables */

  qqp->ind_max_int =qp->nb_int;
  qqp->ind_max_x =qp->n;
  qqp->ind_max_y = qqp->ind_max_x + qp->nb_int*(qp->n) + qp->nb_int*(qp->n -qp->nb_int);
  qqp->ind_max_z =  qqp->ind_max_y + (qp->n)*qp->N;
  qqp->ncol = qqp->ind_max_z + qp->N +1 ;

  int r;
 
  
  r=qqp->ind_dec_x;

  for(i=0;i<r;i++)
    qqp->cont[i]=i;
  
  /* y_ij = sum_k 2^k z_ijk*/
  for(i=0;i<qp->nb_int;i++)
    {
      for (j=0;j<qp->n;j++)
	if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->cont[qqp->ind_dec_x + i*(qp->n) +j] = r;
	    r++;
	  }
    }

  int cpt =0;
  int cpt_bis =0;
  
  /* z_ijk <= u_j_t_ik */
  r= qqp->ind_dec_x + qqp->nb_dec_y;
     
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for(j=0;j<qp->n;j++)
	{
	  for(k=0;k<lim;k++)
	    if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		qqp->cont[qqp->ind_dec_y+ cpt_bis + j*lim + k+cpt] = r;
		r++;
	      }
	}
      cpt=cpt+(qp->n  - qp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qp->nb_int)*lim;
    }
      
  cpt =0;
  cpt_bis =0;
  /* z_ijk <= x_j */
  r= qqp->ind_dec_x + qqp->nb_dec_y+qqp->nb_z_1 ;
     
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for(j=0;j<qp->n;j++)
	{
	  for(k=0;k<lim;k++)
	    if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		qqp->cont[qqp->ind_z_1+ cpt_bis + j*lim + k+cpt] = r;
		r++;
	      }
	}
      cpt=cpt+(qp->n  - qp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qp->nb_int)*lim;
    }


  /* z_ijk >= x_j + u_j_t_ik - u_j */
  r= qqp->ind_dec_x + qqp->nb_dec_y +qqp->nb_z_1 + qqp->nb_z_2;
  cpt =0;
  cpt_bis =0;  
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for(j=0;j<qp->n;j++)
	{
	  for(k=0;k<lim;k++)
	    if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		qqp->cont[qqp->ind_z_2+ cpt_bis + j*lim + k+cpt] = r;
		r++;
	      }
	}
      cpt=cpt+(qp->n  - qp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qp->nb_int)*lim;
    }


  /* z_ijk >= 0 */
  r= qqp->ind_dec_x + qqp->nb_dec_y + qqp->nb_z_1 + qqp->nb_z_2+ qqp->nb_z_3;
  cpt =0;
  cpt_bis =0;   
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for(j=0;j<qp->n;j++)
	{
	  for(k=0;k<lim;k++)
	    if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		qqp->cont[qqp->ind_z_3+ cpt_bis + j*lim + k+cpt] = r;
		r++;
	      }
	}
      cpt=cpt+(qp->n  - qp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qp->nb_int)*lim;
    }

  /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
  r=  qqp->ind_dec_x+ qqp->nb_dec_y+qqp->nb_z_1+qqp->nb_z_2+qqp->nb_z_3+qqp->nb_z_4;
  for(i=0;i<qp->nb_int;i++)
    {
      for (j=i;j<qp->n;j++)
	if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->cont[qqp->ind_z_4 + i* (2*(qp->n) -i +1)/2 + j-i] = r;
	    r++;
	  }
    }

  /*y_ii >= x_i*/ 
  r= qqp->ind_dec_x+ qqp->nb_dec_y+qqp->nb_z_1+qqp->nb_z_2+qqp->nb_z_3+qqp->nb_z_4 + qqp->nb_y_1;
  for(i=0;i<qp->nb_int;i++)
    {
      if(!Zero_BB(qqp->nb_var_y[ij2k(i,i,qp->n)]))
	{
	  qqp->cont[qqp->ind_y_1 + i] = r;
	  r++;
	}
    }
    
  /*y_ij = y_ji*/ 
  r= qqp->ind_dec_x+ qqp->nb_dec_y+qqp->nb_z_1+qqp->nb_z_2+qqp->nb_z_3+qqp->nb_z_4+ qqp->nb_y_1+ qqp->nb_y_2;
  for(i=0;i<qp->nb_int;i++)
    {
      for (j=i+1;j<qp->n;j++)
	if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->cont[qqp->ind_y_2 +i* (2*(qp->n ) -i-1)/2 + j-i-1] = r;
	    r++;
	  }
    }

  /* y_ij <= u_ix_j  for x_i integer and x_j real*/
  r= qqp->ind_dec_x+ qqp->nb_dec_y+qqp->nb_z_1+qqp->nb_z_2+qqp->nb_z_3+qqp->nb_z_4+ qqp->nb_y_1+ qqp->nb_y_2+ qqp->nb_y_3;
  for(i=0;i<qp->nb_int;i++)
    {
      for (j=qp->nb_int;j<qp->n;j++)
	if(!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->cont[qqp->ind_y_3 + i*(qp->n  - qp->nb_int)+ j-qp->nb_int] = r;
	    r++;
	  }
    }

 /* real number of variables +1 for constant term*/	      
  qqp->nb_col =qqp->ind_max_x  + qqp->nb_y + qqp->nb_z + qp->N;
 

  /*nb occurence of each x in constraints involving y or z variables*/
  qqp->x_cont=alloc_vector(qqp->ind_max_x);
  for(i=0;i<qqp->ind_max_x; i++)
    /* initialization to 1 because of the binary expansion constraint for integer x_i*/
    if (i<qqp->nb_int)
      qqp->x_cont[i]=1;
    else
      qqp->x_cont[i]=0;
 

  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      if (!Zero_BB(qqp->nb_var_y[ij2k(i,i,qp->n)]))
	qqp->x_cont[i] = qqp->x_cont[i]+1;
      
      for(j=0;j<qp->n;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	    {
	      qqp->x_cont[i] = qqp->x_cont[i]+1;
	      
	      for(k=0;k<lim;k++)
		qqp->x_cont[j] = qqp->x_cont[j]+2;
	      
	      if (j >= qp->nb_int)
		qqp->x_cont[j] = qqp->x_cont[j]+1;
	    }
	}
    } 

   for(i=qp->nb_int;i<qqp->ind_max_x;i++)
     {
       for(j=0;j<qp->nb_int;j++)
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	    {
	      qqp->x_cont[i] = qqp->x_cont[i]+1;
	    }
     }

  /*table with corresponding variables*/
  qqp->var = alloc_vector(qqp->ncol);

  /*variables x all created*/
  for(i=0;i<qqp->ind_max_x;i++)
    { 
      qqp->var[i]=i;
    }
  
  /* variables y and z*/
  for(i;i<qqp->ncol-qp->N;i++)
    { 
      qqp->var[i]=-1;
    }
  
  /*variables y*/
  r=qqp->ind_max_x;
  
  for(i=0;i<qp->nb_int;i++)
    {
      for (j=0;j<qp->n;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->var[(qp->n)*i + j + qqp->ind_max_x] =r;
	    r++;
	  }
    } 
  
  for(i=qp->nb_int;i<qqp->ind_max_x;i++)
    {
      for (j=0;j<qp->nb_int;j++)
	if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    qqp->var[(qp->n)*qp->nb_int+ (i-qp->nb_int)*qp->nb_int + j +qqp->ind_max_x] =r;
	    r++;
	  }
    }

  /*variables z*/
  cpt=0;
  cpt_bis=0;
  for(i=0;i<qp->nb_int;i++)
    {
      lim=log_base_2(qp->u[i])+1;
      for(j=0;j<qp->n;j++)
	{
	  if (!Zero_BB(qqp->nb_var_y[ij2k(i,j,qp->n)]))
	    for(k=0;k<lim;k++)
	      {
		qqp->var[qqp->ind_max_y +  cpt_bis+j*lim + k+cpt] =r;
		r++;
	      }
	}
      cpt=cpt+(qp->n  - qp->nb_int)*lim;
      cpt_bis=cpt_bis+ (qp->nb_int)*lim;
    }

  /*variables t all created*/
  for(i=qqp->ind_max_z;i<qqp->ncol;i++)
    { 
      qqp->var[i]=r;
      r++;
    }
  
}

/******************************************************************************/
/******************* Structures for branch and bound **************************/
/******************************************************************************/


/******************************************************************************/
/******************************* Create C_MIQCP ********************************/
/******************************************************************************/


C_MIQCP new_c_miqcp()
{
  return (C_MIQCP) malloc(sizeof(struct c_miqcp));
}



void reduce_c_miqcp_obbt( int num_ineg){
  int new_p,new_pq;;
  double * zd;
  double * ze;
  double * zdq;
  double * zeq;
  int i,j,k;
  
  if (num_ineg -cqp->m < cqp->p)
    {
      new_p=cqp->p-1;
      if(!Zero_BB(new_p))
	{
	  zd= alloc_matrix_d(new_p,cqp->n);
	  ze =alloc_vector_d(new_p);

	  for (k=0;k<new_p;k++)
	    if (k + cqp->m < num_ineg)
	      for (i=0;i<cqp->n;i++)
		zd[ij2k(k,i,cqp->n)] = cqp->d[ij2k(k,i,cqp->n)];
	    else
	      if (k + cqp->m >= num_ineg)
		for (i=0;i<cqp->n;i++)
		  zd[ij2k(k,i,cqp->n)] = cqp->d[ij2k(k+1,i,cqp->n)];

	  for (k=0;k<new_p;k++)
	    if (k + cqp->m < num_ineg)
	      ze[k] = cqp->e[k];
	    else
	      if (k + cqp->m >= num_ineg)
		ze[k] = cqp->e[k+1];
	  
	  cqp->p = new_p;
	  free(cqp->d);
	  free(cqp->e);
	  cqp->d=copy_matrix_d(zd,new_p,cqp->n);
	  cqp->e=copy_vector_d(ze,new_p);
	  free(zd);
	  free(ze);
	}
    }
  else
    if (num_ineg -cqp->m-cqp->p-cqp->mq < cqp->pq)
    {
      new_pq=cqp->pq-1;
      if(!Zero_BB(new_pq))
	{
	  zdq= alloc_matrix_3d(new_pq,cqp->n+1,cqp->n+1);
	  zeq =alloc_vector_d(new_pq);

	  for (k=0;k<new_pq;k++)
	    if (k + cqp->m +cqp->p +cqp->mq < num_ineg)
	      for (i=0;i<cqp->n+1;i++)
		for (j=0;j<cqp->n+1;j++)
		  zdq[ijk2l(k,i,j,cqp->n+1,cqp->n+1)] = cqp->dq[ijk2l(k,i,j,cqp->n+1,cqp->n+1)];
	    else
	      if (k + cqp->m +cqp->p +cqp->mq >= num_ineg)
		for (i=0;i<cqp->n;i++)
		  for (j=0;j<cqp->n+1;j++)
		    zdq[ijk2l(k,i,j,cqp->n+1,cqp->n+1)] = cqp->dq[ijk2l(k+1,i,j,cqp->n+1,cqp->n+1)];

	  for (k=0;k<new_pq;k++)
	    if (k + cqp->m +cqp->p +cqp->mq  < num_ineg)
	      zeq[k] = cqp->eq[k];
	    else
	      if (k + cqp->m +cqp->p +cqp->mq  >= num_ineg)
		zeq[k] = cqp->eq[k+1];

	  cqp->pq = new_pq;
	  free(cqp->dq);
	  free(cqp->eq);
	  cqp->dq=copy_matrix_3d(zdq,new_pq,cqp->n+1,cqp->n+1);
	  cqp->eq=copy_vector_d(zeq,new_pq);
	  free(zdq);
	  free(zeq);
	}
     
    }

 
  
}



 



void create_c_miqcp(int flag)
{
  int i,j,l;
  
  cqp = new_c_miqcp();

  double * Q_S = alloc_matrix_d(qp->n,qp->n);
  double * sum_r_Q_r = alloc_matrix_d(qp->n,qp->n);


  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      sum_r_Q_r[ij2k(i,j,qp->n)]=0;
  
  if (flag ==2)
    {
      for(l=0;l<qp->mq;l++)
	for(i=0;i<qp->n;i++)
	  for(j=0;j<qp->n;j++)
	    if (!Zero_BB(qp->aq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)]))
	      sum_r_Q_r[ij2k(i,j,qp->n)]= sum_r_Q_r[ij2k(i,j,qp->n)] + psdp->alphaq[l]*qp->aq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)];
      
  
      for(l=0;l<qp->pq;l++)
	for(i=0;i<qp->n;i++)
	  for(j=0;j<qp->n;j++)
	    if (!Zero_BB(qp->dq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)]))
	      sum_r_Q_r[ij2k(i,j,qp->n)]= sum_r_Q_r[ij2k(i,j,qp->n)] + psdp->alphabisq[l]*qp->dq[ijk2l(l,i+1,j+1,qp->n+1,qp->n+1)];
    }

  /*Calcul de Q - S */
  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      {
	Q_S[ij2k(i,j,qp->n)] = - psdp->beta[ij2k(i,j,qp->n)] -  sum_r_Q_r[ij2k(i,j,qp->n)];
      }


  cqp->N=qp->N;
  cqp->n=qp->n;
  cqp->nb_int=qp->nb_int;
  cqp->m=qp->m;
  cqp->p=qp->p;
  cqp->mq=qp->mq;
  cqp->pq=qp->pq;
  cqp->cons = qp->cons;
  cqp->u = copy_vector_d(qp->u,qp->n + (qp->n+1)*qp->n/2);
  cqp->l = copy_vector_d(qp->l,qp->n + (qp->n+1)*qp->n/2);
   
   
  cqp->q = copy_matrix_d(qp->q, qp->n,qp->n);
  cqp->c = copy_vector_d(qp->c,qp->n);

  if (qp->local_sol == NULL)
    {
      cqp->local_sol=alloc_vector_d(cqp->n);
      cqp->local_sol_adm=MAX_SOL_BB;
    }
  else
    {
      cqp->local_sol=copy_vector_d(qp->local_sol,qp->n);
      cqp->local_sol_adm= qp->sol_adm;
    }
      
  if(!Zero_BB(qp->m))
    {
      cqp->a = copy_matrix_d(qp->a, qp->m,qp->n);
      cqp->b = copy_vector_d(qp->b,qp->m);
    }  
  else
    {
      cqp->a=alloc_matrix_d(1,1);
      cqp->b = copy_vector_d(qp->b,qp->m);
    }
  
  if(!Zero_BB(qp->p))
    {
      cqp->d= copy_matrix_d(qp->d, qp->p,qp->n);
      cqp->e = copy_vector_d(qp->e,qp->p);
    }
  else
    {
      cqp->d=alloc_matrix_d(1,1);
      cqp->e = copy_vector_d(qp->e,qp->p);
    }
  
  if(!Zero_BB(qp->mq))
    {
      cqp->aq = copy_matrix_3d(qp->aq, qp->mq,qp->n+1,qp->n+1);
      cqp->bq = copy_vector_d(qp->bq,qp->mq);
    }  
  else
    {
      cqp->aq=alloc_matrix_d(1,1);
      cqp->bq = copy_vector_d(qp->bq,qp->mq);
    }
  
  if(!Zero_BB(qp->pq))
    {
      cqp->dq= copy_matrix_3d(qp->dq, qp->pq,qp->n+1,qp->n+1);
      cqp->eq = copy_vector_d(qp->eq,qp->pq);
    }
  else
    {
      cqp->dq=alloc_matrix_d(1,1);
      cqp->eq = copy_vector_d(qp->eq,qp->pq);
    }
  
  
   if(!Zero_BB(qp->mq))
     cqp->alphaq = copy_vector_d(psdp->alphaq,qp->mq);
   else
     cqp->alphaq = alloc_vector_d(1);
     
   if(!Zero_BB(qp->pq))
     cqp->alphabisq =copy_vector_d(psdp->alphabisq,qp->pq);
   else
     cqp->alphabisq = alloc_vector_d(1);

   cqp->beta = copy_matrix_d(psdp->beta,qp->n,qp->n);

   int cont = (cqp->n+1)*(cqp->n)/2;
  
   int cont_bis = (cqp->n-1)*(cqp->n)/2;
  
  /*Computation of the initial number of constraints */
  /* Ax=b */
  cqp->ind_eg = qp->m;
  /* Dx = e */
  cqp->ind_ineg = cqp->ind_eg + qp->p;
  /*x^TAqx = bq*/
  cqp->ind_eg_q = cqp->ind_ineg + qp->mq;
  /*x^TDqx = eq*/
  cqp->ind_ineg_q = cqp->ind_eg_q + qp->pq;
  /* y_ij <= u_jx_i + l_ix_j - u_jl_j*/
  if (flag ==3)
    cqp->ind_y_1= cont_bis;
  else
    if (flag ==1)
      cqp->ind_y_1=cont;
    else
      cqp->ind_y_1=cqp->ind_ineg_q + cont;

  /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
  cqp->ind_y_2=cqp->ind_y_1 + cont_bis;

  /* y_ij >= u_ix_j + u_jx_i - u_iu_j */
  if (flag ==3)
    cqp->ind_y_3=cqp->ind_y_2 +cont_bis;
  else
    cqp->ind_y_3=cqp->ind_y_2 +cont;

  /* y_ij >= l_jx_i + l_ix_j - l_il_j */
  if (flag ==3)
    cqp->ind_y_4= cqp->ind_y_3 + cont_bis;
  else
    cqp->ind_y_4= cqp->ind_y_3 + cont;

  /* y_ii >= x_i  or y_ii <= x_i or y_ii = x_i*/
  cqp->nrow= cqp->ind_y_4 + nb_var_int + nb_var_01 + nb_cont_01 + nb_cont_sup_1;
  
 
  int nb_init_cont;
  nb_init_cont = cqp->ind_ineg_q ;

  /* initialisation of the new number of constraints*/
  /* Ax=b */
  cqp->nb_eg =qp->m;
  /* Dx = e */
  cqp->nb_ineg = qp->p;
  /*x^TAqx = bq*/
  cqp->nb_eg_q = qp->mq;
  /*x^TDqx = eq*/
  cqp->nb_ineg_q = qp->pq;
  /* y_ij <= u_jx_i + l_ix_j - u_jl_j*/
  cqp->nb_y_1=0;
  /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
  cqp->nb_y_2=0;
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j */
  cqp->nb_y_3=0;
  /* y_ij >= l_jx_i + l_ix_j - l_il_j */
  cqp->nb_y_4= 0;
  /* y_ii >= x_i  or y_ii <= x_i if flag ==1 */
  cqp->nb_y_5= 0;
  
  cqp->nb_var_y = alloc_matrix(qp->n,qp->n);
  
  for(i=0;i<qp->n;i++)
    for(j=i;j<qp->n;j++)
      {
	/*initialisation*/
	cqp->nb_var_y[ij2k(i,j,qp->n)] =0;
	if (Positive_BB(Q_S[ij2k(i,j,qp->n)]))
	  {
	    cqp->nb_var_y[ij2k(i,j,qp->n)] =1;
	    cqp->nb_var_y[ij2k(j,i,qp->n)] =1;
	  }
	if (Negative_BB(Q_S[ij2k(i,j,qp->n)]))
	  {
	    cqp->nb_var_y[ij2k(i,j,qp->n)] =-1;
	    cqp->nb_var_y[ij2k(j,i,qp->n)] =-1;
	  }
      }

  if (flag ==2)
    {
      for(l=0;l<cqp->mq;l++)
	for(i=0;i<cqp->n;i++)
	  for(j=0;j<cqp->n;j++)
	    if (!Zero_BB(cqp->aq[ijk2l(l,i+1,j+1,cqp->n+1,cqp->n+1)]))
	      cqp->nb_var_y[ij2k(i,j,cqp->n)] =1;
      
      for(l=0;l<cqp->pq;l++)
	for(i=0;i<cqp->n;i++)
	  for(j=0;j<cqp->n;j++)
	    if (!Zero_BB(cqp->dq[ijk2l(l,i+1,j+1,cqp->n+1,cqp->n+1)]) )
	      cqp->nb_var_y[ij2k(i,j,cqp->n)] =1;
    }
  
  cqp->nb_y = nb_non_zero_matrix_sup_i(cqp->nb_var_y, qp->n, qp->n);
  
  
   


  /*computation of the new number of constraints */
  for(i=0;i<qp->n;i++)
    {
      for (j=i;j<qp->n;j++)
	{
	  if (flag == 2)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) )
		{
		  if (i==j)
		    {
		      cqp->nb_y_1++;
		      cqp->nb_y_3++;
		      cqp->nb_y_4++;
		    }
		  else
		    {
		      cqp->nb_y_1++;
		      cqp->nb_y_2++;
		      
		      cqp->nb_y_3++;
		      cqp->nb_y_4++;
		    }
		  if ((i==j) && ( (i<cqp->nb_int) || (i >= cqp->nb_int && ((cqp->u[i]==1 && cqp->l[i]==0) || cqp->l[i]>=1))))
		    cqp->nb_y_5++;
		}
	    }
	  
	  if ((flag == 1) || (flag==3))
	    {
	      if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) )
		{
		  if ((i==j) && (flag==1))
		    cqp->nb_y_1++;
		  else
		    {
		      if (i!=j)
			{
			  cqp->nb_y_1++;
			  cqp->nb_y_2++;
			}
		    }
		}
	      
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) )
		{
		  if ((i!=j) && (flag ==3))
		    {
		      cqp->nb_y_3++;
		      cqp->nb_y_4++;
		    }
		  if (flag==1)
		    {
		      cqp->nb_y_3++;
		      cqp->nb_y_4++;
		    }
		}
	      
	      if ((i==j)  &&  ( (Negative_BB(cqp->nb_var_y[ij2k(i,i,qp->n)]) && (i >= cqp->nb_int && ((cqp->u[i]==1 && cqp->l[i]==0) || cqp->l[i]>=1))) || (!Zero_BB(cqp->nb_var_y[ij2k(i,i,qp->n)]) && (i < qp->nb_int))))
		cqp->nb_y_5++;
	    }
	}
    }

    if (flag == 1 || flag ==3)
      cqp->nb_row =  cqp->nb_y_1+  cqp->nb_y_2 + cqp->nb_y_3+ cqp->nb_y_4+ cqp->nb_y_5;
    else
      cqp->nb_row = cqp->nb_eg +cqp->nb_ineg + cqp->nb_eg_q  +cqp->nb_ineg_q +   cqp->nb_y_1+  cqp->nb_y_2 + cqp->nb_y_3+ cqp->nb_y_4+ cqp->nb_y_5;
  
  
  
  cqp->cont = alloc_vector(cqp->nrow);
  for(i=0;i<cqp->nrow;i++)
    cqp->cont[i]=-1;

 
  
  int r;
  r=nb_init_cont;
 
  for(i=0;i<r;i++)
    cqp->cont[i]=i;
  
  /* y_ij <= l_ix_j + u_jx_i - l_iu_j*/
  for(i=0;i<qp->n;i++)
    {
      for (j=i;j<qp->n;j++)
	{
	  if (flag==2)
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		cqp->cont[cqp->ind_ineg_q + i* (2*(qp->n ) -i +1)/2 + j-i] =  r;
		r++;
	      }
	  if (flag==1)
	    {
	       if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		 {
		   cqp->cont[ i* (2*(qp->n ) -i +1)/2 + j-i] =  r;
		   r++;
		 }
	    }
	  if (flag==3)
	    {
	      if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) && (j>i))
		 {
		   cqp->cont[ i* (2*(qp->n ) -i -1)/2 + j-i-1] =  r;
		   r++;
		 }
	    }
	}
    }
    

  /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
  r= nb_init_cont+ cqp->nb_y_1;
     
  for(i=0;i<qp->n;i++)
    {
      for (j=i+1;j<qp->n;j++)
	{
	  if (flag==2)
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		cqp->cont[cqp->ind_y_1 + i* (2*(qp->n ) - i -1)/2 + j-i-1]=  r;
		r++;
	      }
	  if ((flag==1) || (flag ==3))
	    {
	      if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{
		  cqp->cont[cqp->ind_y_1 + i* (2*(qp->n ) - i -1)/2 + j-i-1]=  r;
		  r++;
		}
	    }
	}
    }

  /* y_ij >= u_ix_j + u_jx_i - u_iu_j*/ 
  r= nb_init_cont+ cqp->nb_y_1+cqp->nb_y_2;
  
  for(i=0;i<qp->n;i++)
    {
      for (j=i;j<qp->n;j++)
	{
	  if (flag==2)
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		cqp->cont[cqp->ind_y_2 + i* (2*(qp->n ) -i +1)/2 + j-i] = r;
		r++;
	      }
	  if (flag==1)
	    {
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{
		  cqp->cont[cqp->ind_y_2 + i* (2*(qp->n ) -i +1)/2 + j-i] = r;
		  r++;
		}
	    }
	  if (flag==3)
	    {
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) && (j>i))
		{
		  cqp->cont[cqp->ind_y_2 + i* (2*(qp->n ) -i -1)/2 + j-i-1] = r;
		  r++;
		}
	    }
	
	}
    }



  
  /* y_ij >= l_ix_j + l_jx_i - l_il_j*/ 
  r= nb_init_cont+ cqp->nb_y_1+cqp->nb_y_2+cqp->nb_y_3;
    
  for(i=0;i<qp->n;i++)
    {
      for (j=i;j<qp->n;j++)
	{
	  if (flag==2)
	    if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
	      {
		cqp->cont[cqp->ind_y_3 + i* (2*(qp->n ) - i +1)/2 + j-i]= r;
		r++;
	      }
	  if (flag==1)
	    {
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{
		  cqp->cont[cqp->ind_y_3 + i* (2*(qp->n ) - i +1)/2 + j-i]= r;
		  r++;
		}
	    }
	  if (flag==3)
	    {
	      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]) && (j>i))
		{
		  cqp->cont[cqp->ind_y_3 + i* (2*(qp->n ) - i -1)/2 + j-i-1]= r;
		  r++;
		}
	    }
	}
    }
  
      
  /* y_ii >= x_i (integer or cont with l[i]>=1) y_ii <= x_i if cont and between 0 and 1, or y_ii = x_i (case binary)*/ 
  r= nb_init_cont+ cqp->nb_y_1+cqp->nb_y_2+cqp->nb_y_3+cqp->nb_y_4;
      
  for(i=0;i<qp->n;i++)
    {
      if (((flag==2) && ( (i<cqp->nb_int) || (i >= cqp->nb_int && ((cqp->u[i]==1 && cqp->l[i]==0) || cqp->l[i]>=1))) ) || flag == 3)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,i,qp->n)]))
	  {
	    cqp->cont[cqp->ind_y_4 + i]= r;
	    r++;
	  }	 
		 
      if (flag==1  &&  ( (Negative_BB(cqp->nb_var_y[ij2k(i,i,qp->n)]) && (i >= cqp->nb_int)) || (!Zero_BB(cqp->nb_var_y[ij2k(i,i,qp->n)]) && (i < qp->nb_int))))
	{
	  cqp->cont[cqp->ind_y_4 + i]= r;
	  r++;
	}
    }
     

 
  /* number of initial variables */

  cqp->ind_max_int =qp->nb_int;
  cqp->ind_max_x =qp->n ;
  cqp->ncol =cqp->ind_max_x + qp->n*qp->n ;
  

  /* real number of variables*/
  cqp->nb_col =cqp->ind_max_x + cqp->nb_y+1;
  
  /*nb occurence of each x in constraints */
  cqp->x_cont=alloc_vector(cqp->ind_max_x);
  for(i=0;i<cqp->ind_max_x; i++)
    cqp->x_cont[i]=0;
 
  
  for(i=0;i<cqp->ind_max_x;i++)
    {
      for (j=i;j<cqp->ind_max_x;j++)
	{
	  if (flag==2)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{ 
		  if ((i==j) && ((i<cqp->nb_int) || (i >= cqp->nb_int && ((cqp->u[i]==1 && cqp->l[i]==0) || cqp->l[i]>=1))) )
		    {
		      cqp->x_cont[i] = cqp->x_cont[i]+4;
		    }
		  else
		    if (i==j)
		      {
			cqp->x_cont[i] = cqp->x_cont[i]+3;
		      }
		    else
		      {
			cqp->x_cont[i] = cqp->x_cont[i]+4;
			cqp->x_cont[j] = cqp->x_cont[j]+4;
		      }
		}
	    }
	  if (flag ==1)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{
		  if (i==j && i < cqp->nb_int)
		    cqp->x_cont[i] = cqp->x_cont[i]+1;
		  else
		    {
		      if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
			{
			  if (i==j)
			    cqp->x_cont[i] = cqp->x_cont[i]+2;
			  else
			    {
			      cqp->x_cont[i] = cqp->x_cont[i]+2;
			      cqp->x_cont[j] = cqp->x_cont[j]+2;
			    }
			}
		      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
			{ 
			  cqp->x_cont[i] = cqp->x_cont[i]+2;
			  cqp->x_cont[j] = cqp->x_cont[j]+2;
			}
		    }
		}
	    }
	 
	  if (flag ==3)
	    {
	      if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
		{
		  if (i==j)
		    {
		      cqp->x_cont[i] = cqp->x_cont[i]+1;
		    }
		  else
		    {
		      if (Negative_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
			{ 
			  cqp->x_cont[i] = cqp->x_cont[i]+2;
			  cqp->x_cont[j] = cqp->x_cont[j]+2;
			}
		      
		      if (Positive_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
			{ 
			  cqp->x_cont[i] = cqp->x_cont[i]+2;
			  cqp->x_cont[j] = cqp->x_cont[j]+2;
			}
		    }
		}
	    }

	
	}
      
    }

  /*table with corresponding variables*/
  cqp->var = alloc_vector(cqp->ncol);

  for(i=0;i<cqp->ind_max_x;i++)
    { 
      cqp->var[i]=i;
    }
  
  for(i;i<cqp->ncol;i++)
    { 
      cqp->var[i]=-1;
    }

  r=cqp->ind_max_x;
  
  for(i=0;i<qp->n;i++)
    {
      for (j=i;j<qp->n;j++)
	if (!Zero_BB(cqp->nb_var_y[ij2k(i,j,qp->n)]))
	  {
	    cqp->var[(qp->n)*i + j +cqp->ind_max_x] =r;
	    r++;
	  }
	}
  
  
  
}

/******************************************************************************/
/******************************* Create MBAB **********************************/
/******************************************************************************/
MIQCP_BAB create_mbab_liste()
{
  MIQCP_BAB bab = (MIQCP_BAB)malloc(sizeof(struct miqcp_bab));
  return bab;
}

struct miqcp_bab create_mbab()
{
  int i,j,k;
  struct miqcp_bab bab;
  
  bab.n=qp->n;
  bab.p=qp->p;
  bab.pq=qp->pq;
  bab.mq=qp->mq;
  bab.nb_int=qp->nb_int;
  
  bab.u=alloc_vector_d(cqp->nb_col);
  bab.l=alloc_vector_d(cqp->nb_col);

  for(i=0;i<bab.n;i++)
    {
      bab.u[i] = qp->u[i];
      bab.l[i] = qp->l[i];
    }

  
  
  /* Attention cas general le faire pour les 4 produits possible*/
  for(j=0;j<bab.n;j++)
    for(k=j;k<bab.n;k++)
      	{
	  if (!Zero_BB(cqp->nb_var_y[ij2k(j,k,bab.n)]))
	    {
	      if (k==j)
		{
		  if ((Positive_BB(qp->u[j])||Zero_BB(qp->u[j])) && (Positive_BB(qp->l[j]) || Zero_BB(qp->l[j])))
		    {
		      bab.u[i]=qp->u[j]*qp->u[j];
		      bab.l[i]=qp->l[j]*qp->l[j];
		      i++;
		    }
		  else
		    {
		      if ((Positive_BB(qp->u[j])|| Zero_BB(qp->u[j])) && Negative_BB(qp->l[j]))
			{
			  if (qp->u[j]*qp->u[j] >qp->l[j]*qp->l[j]) 
			    bab.u[i]=qp->u[j]*qp->u[j];
			  else
			    bab.u[i]=qp->l[j]*qp->l[j];
			  bab.l[i]=qp->l[j]*qp->u[j];
			  i++;
			}
		      else
			if (Negative_BB(qp->u[j]) && Negative_BB(qp->l[j]))
			  {
			    bab.u[i]=qp->l[j]*qp->l[j];
			    bab.l[i]=qp->u[j]*qp->u[j];
			    i++;
			  }
		    }
		}
	      else
		{
		  if ((Positive_BB(qp->u[j])||Zero_BB(qp->u[j])) && (Positive_BB(qp->u[k]) || Zero_BB(qp->u[k])))
		    {
		      if (Negative_BB(qp->l[j]) && Negative_BB(qp->l[k]))
			{
			  if (qp->u[j]*qp->u[k] > qp->l[j]*qp->l[k])
			    bab.u[i]=qp->u[j]*qp->u[k];
			  else
			    bab.u[i]=qp->l[j]*qp->l[k];
			  if  (qp->l[j]*qp->u[k] > qp->l[k]*qp->u[j])
			    bab.l[i]=qp->l[j]*qp->u[k];
			  else
			    bab.l[i]=qp->l[k]*qp->u[j];
			  i++;
			}
		      
		      else
			{
			  if (Negative_BB(qp->l[j]))
			    {
			      bab.u[i]=qp->u[j]*qp->u[k];
			      bab.l[i]=qp->l[j]*qp->u[k];
			      i++;
			    }
			  else
			    {
			      if (Negative_BB(qp->l[k]))
				{
				  bab.u[i]=qp->u[j]*qp->u[k];
				  bab.l[i]=qp->u[j]*qp->l[k];
				  i++;
				}
			      else
				{
				  bab.u[i]=qp->u[j]*qp->u[k];
				  bab.l[i]=qp->l[j]*qp->l[k];
				  i++;
				}
			    }
			}
		    }
		  else
		    {
		      if (Negative_BB(qp->u[j]) && Negative_BB(qp->u[k]) )
			{
			  bab.u[i]=qp->l[k]*qp->l[j];
			  bab.l[i]=qp->u[k]*qp->u[j];
			  i++;
			}
		      else
			{
			  if (Negative_BB(qp->u[j]))
			    {
			      bab.u[i]=qp->l[k]*qp->l[j];
			      bab.l[i]=qp->u[k]*qp->l[j];
			      i++;
			    }
			  else
			    if (Negative_BB(qp->u[k]))
			      {
				bab.u[i]=qp->l[k]*qp->l[j];
				bab.l[i]=qp->u[j]*qp->l[k];
				i++;
			      }
			  
			}
		    }
		  
		}
	      
	    }
	}
 
	
  bab.i=0;
  bab.j=0;
  return bab;
}



MIQCP_BAB realloc_mbab_liste(struct miqcp_bab bab, int nb_col)
{
  MIQCP_BAB new_bab =  create_mbab_liste();
  new_bab->n=bab.n;
  new_bab->p=bab.p;
  new_bab->pq=bab.pq;
  new_bab->nb_int=bab.nb_int;
    
  new_bab->u=alloc_vector_d(nb_col);
  memcpy(new_bab->u,bab.u,nb_col*sizeof(double));
  
  new_bab->l=alloc_vector_d(nb_col);
  memcpy(new_bab->l,bab.l,nb_col*sizeof(double));
  
  return new_bab;
}

/******************************************************************************/
/******************* Structures for the PSD solver ****************************/
/******************************************************************************/

SDP new_sdp()
{
  return (SDP) malloc(sizeof(struct sdp));
}


/******************************************************************************/
/******************************* Create tab cont*******************************/
/******************************************************************************/

TAB_CONT create_tab_cont()
{
  int i,j;

  TAB_CONT tab;
  
  tab = (int **)calloc(psdp->length,sizeof(int*)); 
    for(i=0;i<psdp->length;i++)
    tab[i]=alloc_vector(3);

  int cont_X_1=psdp->i1;
  /* Constraints  x_ix_j - u_ix_j <= 0*/ 
  for(i=0;i<psdp->nb_int;i++)
    {
      for(j=i+1;j<psdp->n ;j++)
      {
	tab[cont_X_1-1][0]=i;
	tab[cont_X_1-1][1]=j;
	tab[cont_X_1-1][2]=1;
	cont_X_1++;
      }
    }

 /* Constraints x_ix_j - u_jx_i <= 0*/
  int cont_X_2=psdp->i2;
  for(i=0;i<psdp->nb_int;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_2-1][0]=i;
	  tab[cont_X_2-1][1]=j;
	  tab[cont_X_2-1][2]=2;
	  cont_X_2++;
	}
    }

  /* Constraints - x_ix_j + u_jx_i + u_iu_j - u_iu_j <= 0*/
  int cont_X_3=psdp->i3;
  for(i=0;i<psdp->nb_int;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_3-1][0]=i;
	  tab[cont_X_3-1][1]=j;
	  tab[cont_X_3-1][2]=3;
	  cont_X_3++;
	}
    }
  
  /* Constraints - x_ix_j <= 0*/
  int cont_X_4=psdp->i4;
  for(i=0;i<psdp->nb_int;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_4-1][0]=i;
	  tab[cont_X_4-1][1]=j;
	  tab[cont_X_4-1][2]=4;
	  cont_X_4++;
	}
    }
  
  return tab;
 
}

TAB_CONT create_tab_cont_mixed()
{
  int i,j;

  TAB_CONT tab;
  
  tab = (int **)calloc(psdp->length,sizeof(int*)); 
    for(i=0;i<psdp->length;i++)
    tab[i]=alloc_vector(3);

  int cont_X_1=psdp->i1;
  /* Constraints  x_ix_j - u_ix_j <= 0*/ 
  for(i=0;i<psdp->n;i++)
    {
      for(j=i+1;j<psdp->n ;j++)
      {
	tab[cont_X_1-1][0]=i;
	tab[cont_X_1-1][1]=j;
	tab[cont_X_1-1][2]=1;
	cont_X_1++;
      }
    }

 /* Constraints x_ix_j - u_jx_i <= 0*/
  int cont_X_2=psdp->i2;
  for(i=0;i<psdp->n;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_2-1][0]=i;
	  tab[cont_X_2-1][1]=j;
	  tab[cont_X_2-1][2]=2;
	  cont_X_2++;
	}
    }

  /* Constraints - x_ix_j + u_jx_i + u_iu_j - u_iu_j <= 0*/
  int cont_X_3=psdp->i3;
  for(i=0;i<psdp->n;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_3-1][0]=i;
	  tab[cont_X_3-1][1]=j;
	  tab[cont_X_3-1][2]=3;
	  cont_X_3++;
	}
    }
  
  /* Constraints - x_ix_j <= 0*/
  int cont_X_4=psdp->i4;
  for(i=0;i<psdp->n;i++)
    {
      for(j=i+1;j<psdp->n;j++)
	{
	  tab[cont_X_4-1][0]=i;
	  tab[cont_X_4-1][1]=j;
	  tab[cont_X_4-1][2]=4;
	  cont_X_4++;
	}
    }
  
  return tab;
 
}




/******************************************************************************/
/******************************* Create SDP ***********************************/
/******************************************************************************/
void create_sdp()
{
  int i,j,k;
  psdp = new_sdp();

  psdp->n = qp->n;
  psdp->nb_int = qp->nb_int;
  psdp->m = qp->m;
  psdp->p = qp->p;
  psdp->mq = qp->mq;
  psdp->pq = qp->pq;
  psdp->miqp = qp->miqp;
  psdp->cons = -qp->cons;
 
  psdp->q = copy_matrix_d(qp->q,qp->n,qp->n);
  
  psdp->a = copy_matrix_d(qp->a,qp->m,qp->n);
  psdp->d = copy_matrix_d(qp->d,qp->p,qp->n);

  psdp->aq = copy_matrix_3d(qp->aq,qp->mq,qp->n+1,qp->n+1);
  psdp->dq = copy_matrix_3d(qp->dq,qp->pq,qp->n+1,qp->n+1);
  
  psdp->c = copy_vector_d(qp->c,qp->n);
  psdp->u = copy_vector_d(qp->u,qp->n);
  psdp->l = copy_vector_d(qp->l,qp->n);
  psdp->b = copy_vector_d(qp->b,qp->m);
  psdp->e = copy_vector_d(qp->e,qp->p);
  
  psdp->bq = copy_vector_d(qp->bq,qp->mq);
  psdp->eq = copy_vector_d(qp->eq,qp->pq);
  
  psdp->beta=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_1=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_2=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_3=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_4=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_diag=alloc_vector_d(qp->nb_int);
  
  /*initialize beta matrix*/
  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      {
	psdp->beta_1[ij2k(i,j,qp->n)]=0;
	psdp->beta_2[ij2k(i,j,qp->n)]=0;
	psdp->beta_3[ij2k(i,j,qp->n)]=0;
	psdp->beta_4[ij2k(i,j,qp->n)]=0;
      }

  for(i=0;i<qp->nb_int;i++)
    psdp->beta_diag[i]=0;
  
  psdp->alpha=0;
  psdp->alphabis=0;
  psdp->alphaq=alloc_vector_d(qp->mq);
  psdp->alphabisq=alloc_vector_d(qp->pq);
  
  psdp->length=nb_max_cont(psdp,&psdp->i1,&psdp->i2,&psdp->i3,&psdp->i4);
  psdp->nb_cont = FACTOR * psdp->length;
  psdp->nb_max_cont = psdp->nb_cont;
 
  psdp->constraints=alloc_vector(psdp->nb_cont );
  psdp->constraints_old=alloc_vector( psdp->nb_cont );
  for(i=0;i<psdp->nb_cont;i++)
    {
      psdp->constraints[i]=0;
      psdp->constraints_old[i]=0;
    }

  
  psdp->tab=create_tab_cont();

}

void create_sdp_mixed()
{
  int i,j,k;

  psdp = new_sdp();

  psdp->n = qp->n;
  psdp->nb_int = qp->nb_int;
  psdp->m = qp->m;
  psdp->p = qp->p;
  psdp->mq = qp->mq;
  psdp->pq = qp->pq;
  psdp->miqp = qp->miqp;
   psdp->cons = -qp->cons;
 
  psdp->q = copy_matrix_d(qp->q,qp->n,qp->n);
  
  psdp->a = copy_matrix_d(qp->a,qp->m,qp->n);
  psdp->d = copy_matrix_d(qp->d,qp->p,qp->n);

  psdp->aq = copy_matrix_3d(qp->aq,qp->mq,qp->n+1,qp->n+1);
  psdp->dq = copy_matrix_3d(qp->dq,qp->pq,qp->n+1,qp->n+1);
  
  psdp->c = copy_vector_d(qp->c,qp->n);
  psdp->u = copy_vector_d(qp->u,qp->n);
  psdp->l = copy_vector_d(qp->l,qp->n);
  psdp->b = copy_vector_d(qp->b,qp->m);
  psdp->e = copy_vector_d(qp->e,qp->p);
  
  psdp->bq = copy_vector_d(qp->bq,qp->mq);
  psdp->eq = copy_vector_d(qp->eq,qp->pq);
  
  psdp->beta=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_1=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_2=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_3=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_4=alloc_matrix_d(qp->n,qp->n);
  psdp->beta_diag=alloc_vector_d(qp->n);
  
  /*initialize beta matrix*/
  for(i=0;i<qp->n;i++)
    for(j=0;j<qp->n;j++)
      {
	psdp->beta_1[ij2k(i,j,qp->n)]=0;
	psdp->beta_2[ij2k(i,j,qp->n)]=0;
	psdp->beta_3[ij2k(i,j,qp->n)]=0;
	psdp->beta_4[ij2k(i,j,qp->n)]=0;
      }

  for(i=0;i<qp->n;i++)
    psdp->beta_diag[i]=0;
  
  psdp->alpha=0;
  psdp->alphabis=0;
  psdp->alphaq=alloc_vector_d(qp->mq);
  psdp->alphabisq=alloc_vector_d(qp->pq);
  
  psdp->length=nb_max_cont_mixed(psdp,&psdp->i1,&psdp->i2,&psdp->i3,&psdp->i4);
  psdp->nb_cont = FACTOR * psdp->length;
  psdp->nb_max_cont = psdp->nb_cont;
 
  psdp->constraints=alloc_vector(psdp->nb_cont );
  psdp->constraints_old=alloc_vector( psdp->nb_cont );
  for(i=0;i<psdp->nb_cont;i++)
    {
      psdp->constraints[i]=0;
      psdp->constraints_old[i]=0;
    }
 
  
  psdp->tab=create_tab_cont_mixed();

    
}




/******************************************************************************/
/******************************* Free Memory **********************************/
/******************************************************************************/

void free_miqcp()
{
  if(Positive_BB(qp->n))
    {
      free_vector_d(qp->q);
      free_vector_d(qp->c);
      free_vector_d(qp->u);
      free_vector_d(qp->l);
      free_vector(qp->lg);
      free_vector(qp->lgc);
    }
  
  if(Positive_BB(qp->m))
    {
      free_vector_d(qp->a);
      free_vector_d(qp->b);
    }

  if(Positive_BB(qp->p))
    {
      free_vector_d(qp->d);
      free_vector_d(qp->e);
    }
  
  if(Positive_BB(qp->mq))
    {
      free_vector_d(qp->aq);
      free_vector_d(qp->bq);
    }
   
   if(Positive_BB(qp->pq))
    {
      free_vector_d(qp->dq);
      free_vector_d(qp->eq);
    }
   
  free(qp);
}

void free_sdp()
{
  if(Positive_BB(psdp->n))
    {
      free_vector_d(psdp->q);
      free_vector_d(psdp->c);
      free_vector_d(psdp->u);
      free_vector_d(psdp->beta_1);
      free_vector_d(psdp->beta_2);
      free_vector_d(psdp->beta_3);
      free_vector_d(psdp->beta_4);
      free_vector(psdp->constraints);
    }
  
  if(Positive_BB(psdp->m))
    {
      free_vector_d(psdp->a);
      free_vector_d(psdp->b);
    }
  
  if(Positive_BB(psdp->p))
    {
      free_vector_d(psdp->d);
      free_vector_d(psdp->e);
    }

  if(Positive_BB(psdp->mq))
    {
      free_vector_d(psdp->aq);
      free_vector_d(psdp->bq);
    }
 
  if(Positive_BB(psdp->pq))
    {
      free_vector_d(psdp->dq);
      free_vector_d(psdp->eq);
    }
  
    free(psdp);
}


void free_c_miqcp()
{
  if(cqp->n>0)
    {
      free_vector_d(cqp->q);
      free_vector_d(cqp->new_q);
      free_vector_d(cqp->c);
      free_vector_d(cqp->u);
      free_vector_d(cqp->l);
    }
  
  if(Positive_BB(cqp->m))
    {
      free_vector_d(cqp->a);
      free_vector_d(cqp->a_t_a);
      free_vector_d(cqp->b);
    }
  
  if(Positive_BB(cqp->p))
    {
      free_vector_d(cqp->d);
      free_vector_d(cqp->d_t_d);
      free_vector_d(cqp->e);
    }
  
  if(Positive_BB(cqp->mq))
    {
      free_vector_d(cqp->aq);
      free_vector_d(cqp->bq);
    }
  
   if(Positive_BB(cqp->pq))
    {
      free_vector_d(cqp->dq);
      free_vector_d(cqp->eq);
    }

    
  free_vector_d(cqp->beta);
  free(cqp);
}

void free_tab_cont(TAB_CONT tab, int n)
{
  int i;
  for (i=0;i<n;i++)
    free(tab[i]);
}
