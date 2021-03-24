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

#ifndef UTILITIES_H
#define UTILITIES_H



#include"liste.h"
#include"quad_prog.h"
//#include "parameters.h"

extern double EPS_BETA;
extern double EPS_LM;
extern double EPS_VIOL_CB;


/* Macros */
#define Round(a) (a>=0?(int) (a+0.5) :(int) (a-0.5) ) 
#define Positive_BB(a) (a> EPS_BETA) 
#define Negative_BB(a) (a<- EPS_BETA)
#define Zero_BB(a) (!Positive_BB(a) && !Negative_BB(a))
#define Positive_VIOL_CB(a) (a> EPS_VIOL_CB) 
#define Negative_VIOL_CB(a) (a< - EPS_VIOL_CB)
#define Zero_VIOL_CB(a) (!Positive_VIOL_CB(a) && !Negative_VIOL_CB(a))
#define Positive_LM(a) (a> EPS_LM) 
#define Negative_LM(a) (a<- EPS_LM) 
#define Zero_LM(a) (!Positive_LM(a) && !Negative_LM(a))

/*indices from an ** to an *, and reciprocally*/
#define ij2k(i,j,n) ((i)*(n)+(j))
#define k2i(k,n) ((int)(k)/(n))
#define k2j(k,n) ((k)%(n))

/*indices from an *** to an *, and reciprocally*/
#define ijk2l(i,j,k,n,p) ((i)*(n)*(p)+(j)*(p) +(k))


/*  print tables */
void print_vec_d(double * tab, int n);
void print_mat_d(double * tab, int n, int m);
void print_vec(int * tab, int n);
void print_mat(int * tab, int n, int m);
void print_tab(char * tab, int n);
void print_mat_3d(double * tab, int n, int m,int l);
/* allocate memory */
int* alloc_vector (int dimension);
int * alloc_matrix (int ligne, int colonne);
double* alloc_vector_d (int dimension);
double * alloc_matrix_d (int ligne, int colonne);
double * alloc_matrix_3d (int n, int m, int p);
char* alloc_string (int dimension);
char* alloc_tab (int dim1,int dim2);
MIQCP_BAB realloc_mbab_liste(struct miqcp_bab bab, int nb_col);
void free_string(char * chaine);
void free_vector(int * v);
void free_vector_d(double * v);
void free_and_null (char **ptr);

/* basic operations on vectors*/
double v_abs(double a);
void v_abs_ref(double *a);
int sum_vector(int * a, int j,int n);
double sum_vector_d(double * a, int j, int n);
int maximum(int * vector,int n);
double maximum_d(double * vector,int n);
double maximum_d_abs(double * vector,int n);
double minimum_d(double * vector,int n);
int minimum(int * vector,int n);
int * copy_vector(int *v, int n);
double * copy_vector_d(double *v, int n);
double * copy_matrix_d(double *matrix, int m,int n);
double * copy_matrix_3d(double *matrix, int m,int n,int p);
char * copy_string(char *s);
void copy_line(double * B, double * v, int i,int t);
int nb_non_zero_matrix_sym(double * a, int n);
int nb_non_zero_matrix(double * a, int n, int m);
int nb_non_zero_matrix_i(int * a, int n, int m);
int nb_non_zero_matrix_sup_i(int *a, int n, int m);
int nb_non_zero_matrix_without_diag(double *a, int n, int m);
int nb_non_zero_matrix_without_diag_i(int *a, int n, int m);
int nb_non_zero_vector(double * v, int n);
int nb_non_zero_vector_int(int * v, int n);
int nb_non_zero_matrix_3d(double * a, int m, int n, int p);
double * add_two_vector(double * v1, int * v2,int n);
double * opposite_matrix_3d(double * M, int m, int n, int p);

/*management of bounds*/
int log_base_2 (int x);
int vector_length_base_2(double * u, int n);
int * create_u_base2_cumulated(double * u, int n);
int * create_u_base2(double * u, int n);

/*change an int into a char*/
char * intochar(int val);
double random_double(double min, double max);

/* basic operations on matrices*/
double * multiplication_matrices( double * a, double * b, int n);
double * multiplie_d(double * M, double * P, int n,int m);
double * transpose_d(double * M, int n, int m);
double * matrice_par_scal(double * M, int n, int m, double scal);
double * addition_mat_d(double* M,double * P,int n);
double * carre_cont(double * M, int m, int n);
double * vect_c(double* M,double * v, int m, int n);
double * vect_par_scal(double * v, int n, double scal);
double * addition_vect(double * v,double * w, int n);
double constante(double * v, int n,int scal);
int is_symetric(double * a, int n);
int is_symetric_3d(double * a,int k, int n);
double * transpose_matrix( double * a, int n);
double * multiplication_matrices( double * a, double * b, int n);
double * transpose_matrix_mn( double * a, int m, int n);

/* utilities for cplex*/
double sum_r_ari_b_r(double * a, double * b,int m,int n, int i);
double sum_r_ari_arj(double * a, int m,int n, int i, int j);
double sum_r_b_r(double * b,int m);
double sum_r_alphaq_bq_r(double * b,double* alphaq,int m);

/* utilities for the sdp solver*/
void dualized_objective_function(SDP psdp, double ** q_beta, double ** c_beta, double * l);
void evaluate_objective(double * objective_value,SDP psdp, double ** q_beta, double ** c_beta, double *l_beta);
void compute_subgradient(double ** subgrad,SDP psdp);
void compute_new_subgradient(double ** check_violated,SDP psdp,int new_length, int ** variable_indices,double **X);
void purge_constraints_not_violated(SDP psdp, double ** slacks, int * n_assign, int ** assign_new_from_old);
void add_constraints_violated(SDP psdp, double ** check_violated,int * n_append);
void initialize_constraints_violated(SDP psdp, double ** check_violated, int * new_length);
int check_constraints_violated(SDP psdp, double ** check_violated, double ** x);
void update_beta_value(SDP psdp, double * zbeta);
int nb_max_cont(SDP psdp, int * ind_X1,int * ind_X2,int *ind_X3,int*ind_X4);
int nb_cont_Xij(int n, int nb_int, int p);
int nb_max_cont_cont(SDP psdp, int * ind_X1,int * ind_X2,int *ind_X3,int*ind_X4);
int nb_cont_Xij_cont(int n);
int nb_max_cont_mixed(SDP psdp, int * ind_X1,int * ind_X2,int *ind_X3,int*ind_X4);

/*******************************  utilities for sbb **********************************************/


void refine_bounds(struct miqcp_bab bab);

int feasible_for_constraints(double * x_y);

int select_i_j_miqcrq(double * x_y, struct miqcp_bab bab , int *i, int *j);

struct miqcp_bab copy_bab(struct miqcp_bab bab);

struct miqcp_bab update_variables_int_branch_1(struct miqcp_bab bab, double * x_y, int *i,int *j);
struct miqcp_bab update_variables_int_branch_2(struct miqcp_bab bab, double * x_y, int *i,int *j);


struct miqcp_bab update_variables_cont_branch_1(struct miqcp_bab bab, double * x_y, int *i,int *j);
struct miqcp_bab update_variables_cont_branch_2(struct miqcp_bab bab, double * x_y, int *i,int *j);



int test_bound_miqcp(struct miqcp_bab bab, int nb_col);
int is_integer(double a);
double check_value(double a);
double sum_ij_qi_qj_x_ij(double * q, double * x, int n);
double sum_i_ci_x_i(double * c, double * x, int n);

/* Utilities for changing bounds */ 
void change_bounds();
void change_bounds_back();
void change_bounds_gen_int();
void change_bounds_back_gen_int();

/*************** Computation of Hessian and eigen values  ****************/

void compute_new_q_boxqp();
void compute_new_q_miqcr();
void compute_new_lambda_min_miqcr();
void compute_new_q_sbb_miqcr();
void compute_new_lambda_min_sbb_miqcr();

#endif
