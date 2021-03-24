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

#ifndef QUAD_PROG_H
#define QUAD_PROG_H

struct miqcp
{
  int n;
  int nb_int;
  int m;
  int p;
  int mq;
  int pq;
  int N;
  /* flag_const = 1 if there is a constant term to consider, cons is the value of the constant term*/
  double cons;
  double flag_const;
  /* miqp = 1 if there is continuous variables*/
  int miqp;
  /*flag_alg = 1 => algorithms s_miqcrq or p_miqcrq or s_miqcrq2 or n_miqcrq */
  int flag_alg;
  double *q;
  
  double *a;
  double *d;
  double *c;
  double *u;
  double *l;
  double *init_l;
  double *init_u;
  double *b;
  double * e;
  double * aq;
  double * dq;
  double *bq;
  double *eq;
  int *lg;
  int *lgc;
  double * sol;
  double sol_adm;
  double * local_sol;
 
};

typedef struct miqcp * MIQCP;

struct q_miqcp
{
  int n;
  int nb_int;
  int m;
  int p;
  int mq;
  int pq;
  int N;
  int miqp;
  double cons;
  double *q;
  double * new_q;
  double *a;
  double * a_t_a;
  double *d;
  double * d_t_d;
  double *c;
  double *u;
  double *l;
  double *b;
  double * e;
  double * aq;
  double * sum_r_aq_alphaq;
  double * dq;
  double * sum_r_dq_alphabisq;
  double *bq;
  double * sum_r_bq;
  double *eq;
  double * sum_r_eq;
  int *lg;
  int *lgc;
  double alpha;
  double alphabis;
  double *alphaq;
  double *alphabisq;
  double * beta;
  double * beta_1;
  double * beta_2;
  double * beta_3;
  double * beta_4;
  double * beta_c;
  double * delta;
 
  double new_lambda_min;
  
  int * cont;
  int * x_cont;
  int * nb_var_y;
  int * var;
  double * sol;
  double sol_adm;

  int nb_init_var;
  
  /*number of constraints of each type after reduction of the problem*/
  /* Ax=b */
  int nb_eg;
  /* Dx = e */
  int nb_ineg;
   /*x^TAqx = bq*/
  int nb_eg_q;
  /*x^TDqx = eq*/
  int nb_ineg_q;
  /*Dq,y + Dq_0 x + sq = eq*/
  int nb_ineg_quad;
  
  /* y_ii = x_i */
  int nb_y_1;
  /* y_ij <= x_i i <j */
   int nb_y_2;
  /* y_ij <= x_j i <j */
  int nb_y_3;
  /* y_ij >= x_i + x_j -1 i <j */
  int nb_y_4;
  /* y_ij >= 0 i <j */
  
  int nb_row;
  int nb_cut;
  int nb_cont;
  
  /* for integer case */
  /* x = sum_k 2^k t_ik*/
  int nb_dec_x;
  /* y_ij = sum_k 2^k z_ijk*/
  int nb_dec_y;
  /*z_ijk <= u_jt_ik*/
  int nb_z_1;
  /*z_ijk <= x_j*/
  int nb_z_2;
  /*z_ijk >= x_j - u_j(1 - t_ik)*/
  int nb_z_3;
  /*z_ijk >= 0*/
  int nb_z_4;


  
  /* new number of variables in  the problem*/
  int nb_z;
  int nb_y;
  int nb_col;
  int diag;

  /*number of constraints of each type initialy in the problem*/
   /* Ax=b */
  int ind_eg;
  /* Dx = e */
  int ind_ineg;
   /*x^TAqx = bq*/
  int ind_eg_q;
  /*x^TDqx = eq*/
  int ind_ineg_q;
  /* y_ii = x_i */
  int ind_y_1;
  /* y_ij <= x_i i<j*/
  int ind_y_2;
  /* y_ij <= x_j i<j*/
  int ind_y_3;
  /* y_ij >= x_i + x_j -1 */
  int ind_y_4;
  /* y_ij >= 0 */
  
  int nrow;

  /* for the general integer case*/
  /* x = sum_k 2^k t_ik*/
  int ind_dec_x;
  /* y_ij = sum_k 2^k z_ijk*/
  int ind_dec_y;
  /*z_ijk <= u_jt_ik*/
  int ind_z_1;
  /*z_ijk <= x_j*/
  int ind_z_2;
  /*z_ijk >= x_j - u_j(1 - t_ik)*/
  int ind_z_3;
  /*z_ijk >= 0*/
  int ind_z_4;
 



  /* number of variables initialy in  the problem*/
  int ind_max_int;;
  int ind_max_x;
  int ind_max_y;
  int ind_max_z;
  int ncol;

};

typedef struct q_miqcp * Q_MIQCP;

/* structures for branch and bound */
struct miqcp_bab
{
  int n;
  int p;
  int pq;
  int mq;
  int nb_int;
  double * u;
  double * l;
  int i;
  int j;
};

typedef struct miqcp_bab * MIQCP_BAB;

struct c_miqcp
{
  int n;
  int nb_int;
  int m;
  int p;
  int mq;
  int pq;
  int N;
  int miqp;
  double cons;
  double local_sol_adm;
  double * local_sol;
  double *q;
  double * new_q;
  double *a;
  double * a_t_a;
  double *d;
  double * d_t_d;
  double *c;
  double *l;
  double *u;
  double *b;
  double * e;
  double * aq;
  double * sum_r_aq_alphaq;
  double * dq;
  double * sum_r_dq_alphabisq;
  double *bq;
  double * sum_r_bq;
  double *eq;
  double * sum_r_eq;
   double alpha;
  double alphabis;
  double *alphaq;
  double *alphabisq;
  double * beta;
  double * beta_1;
  double * beta_2;
  double * beta_3;
  double * beta_4;
  double * beta_c;
  double new_lambda_min;
 
 
  int * cont;
  int * x_cont;

  int * var;
  int * nb_var_y;
  /*number of constraints of each type after reduction of the problem*/
  /* Ax=b */
  int nb_eg;
  /* Dx = e */
  int nb_ineg;
   /*x^TAqx = bq*/
  int nb_eg_q;
  /*x^TDqx = eq*/
  int nb_ineg_q;
  /* y_ij <= u_jx_i + l_ix_j - u_jl_j*/
  int nb_y_1;
  /* y_ij <= u_ix_j + l_jx_i - u_il_j*/
  int nb_y_2;
  /* y_ij >= u_ix_j + u_jx_i - u_iu_j */
  int nb_y_3;
  /* y_ij >= l_jx_i + l_ix_j - l_il_j */
  int nb_y_4;
  /* y_ii >= x_i */
  int nb_y_5;
  int nb_row;
  int nb_cont;
  
  
  /* new number of variables in  the problem*/
  int nb_y;
  int nb_col;
  

  /*number of constraints of each type initialy in the problem*/
   /* Ax=b */
  int ind_eg;
  /* Dx = e */
  int ind_ineg;
   /*x^TAqx = bq*/
  int ind_eg_q;
  /*x^TDqx = eq*/
  int ind_ineg_q;
  /* y_ii = x_i */
  int ind_y_1;
  /* y_ij <= x_i i<j*/
  int ind_y_2;
  /* y_ij <= x_j i<j*/
  int ind_y_3;
  /* y_ij >= x_i + x_j -1 */
  int ind_y_4;
  /* y_ij >= 0 */
  int nrow;

  
 

  /* number of variables initialy in  the problem*/
  int ind_max_int;;
  int ind_max_x;
  int ind_max_y;
  int ncol;

};


typedef struct c_miqcp * C_MIQCP;

/*  structure for the semi-definite program */
typedef int ** TAB_CONT;

struct sdp
{
  int n;
  int nb_int;
  int m;
  int p;
  int mq;
  int pq;
  int miqp;
  int length;
  int nb_cont;
  int nb_max_cont;
  int i1;
  int i2;
  int i3;
  int i4;
 
  int nb_init_var;
  double cons;
  
  double *q;
  double *a;
  double *d;
  double *c;
  double *u;
  double *l;
  double *b;
  double * e;
  double * aq;
  double * dq;
  double *bq;
  double *eq;
  double * x;
  double * beta;
  double * beta_1;
  double * beta_2;
  double * beta_3;
  double * beta_4;
  double * beta_5;
  double * beta_diag;
 
  double alpha;
  double alphabis;
  double *alphaq;
  double *alphabisq;
  int * constraints;
  int * constraints_old;
  TAB_CONT tab;
  char * solver_sdp;
 
};

typedef struct sdp * SDP;



MIQCP new_miqcp();
C_MIQCP new_c_miqcp();
Q_MIQCP new_q_miqcp();

void create_c_miqcp(int flag  );
void create_q_miqcp();
void create_sdp();
void create_sdp_mixed();
struct miqcp_bab create_mbab();


void update_qp();
void reduce_c_miqcp_obbt( int num_ineg);

void free_miqcp();
void free_c_miqcp(); 
void free_sdp(); 
void free_tab_cont(TAB_CONT tab, int n);

#endif
