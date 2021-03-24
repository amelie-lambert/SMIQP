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
#include<time.h>
#include<math.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>
#include<math.h>


#include"liste.h"
#include"quicksort.h"
#include"quad_prog.h"
#include"utilities.h"

extern MIQCP qp;
extern Q_MIQCP qqp;
extern C_MIQCP cqp;  
extern Liste liste;

extern double ub_bb;
extern double lb_bb;
extern double * best_x;
extern long nodecount;
extern double sol_sdp;
extern int father;

extern double EPS;
extern double EPS_BB;
extern double EPS_BRANCH;
extern double EPS_INT;
extern int NB_MAX_CONT;
extern double GAMMA;



extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info); 
/*****************************************************************************/
/***************************** print tables **********************************/
/*****************************************************************************/

void print_tab(char * tab, int n)
{
  int i;
  for (i=0;i<n*n;i++)
    {
      printf("%c ",tab[i]);
      if (i %n ==0)
	printf("\n");
    }
}

void print_vec(int * tab, int n)
{
  int i;
  for (i=0;i<n;i++)
    printf("%d ",tab[i]);
}


void print_vec_d(double * tab, int n)
{
  int i;
  for (i=0;i<n;i++)
    printf("%.12lf ",tab[i]);
}

void print_mat_d(double * tab, int n, int m)
{
  int i;
  for (i=0;i<m*n;i++)
    {
      printf("%.7lf ",tab[i]);
      if ((i+1) %n ==0)
	printf("\n");
    }
}

void print_mat_3d(double * tab, int n, int m,int l)
{
  int i,j,k;
  
  for (i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
      {
	for(k=0;k<l;k++)
	  printf("%.7lf ",tab[ijk2l(i,j,k,m,l)]);
	printf("\n");
      }
      printf("\n");
    }
}
void print_mat(int * tab, int n, int m)
{
  int i;
  for (i=0;i<m*n;i++)
    {
      printf("%d ",tab[i]);
      if ((i+1) %n ==0)
	printf("\n");
    }
}
/*****************************************************************************/
/*************************** allocate  memory ********************************/
/*****************************************************************************/

int * alloc_vector (int dimension)
{
  int* vector = (int *)calloc(dimension,sizeof(int));
  if (vector == NULL)
    exit(EXIT_FAILURE);
  return vector;
}

int * alloc_matrix (int line, int column)
{
  int i;
  int* matrix = (int *)calloc(line*column,sizeof(int));
  if (matrix == NULL)
	exit(EXIT_FAILURE);
  return matrix;
}

double* alloc_vector_d (int dimension)
{
  double* vector = (double*) calloc(dimension,sizeof(double));
  if (vector == NULL)
    exit(EXIT_FAILURE);
  return vector;
}

double * alloc_matrix_d (int line, int column)
{
  int i;
  double* matrix = (double*) calloc(line*column,sizeof(double));
    if (matrix == NULL)
	exit(EXIT_FAILURE);
    return matrix;
}

double * alloc_matrix_3d (int n, int m, int p)
{
  int i,j;
  double* matrix = (double *)calloc(n*m*p,sizeof(double));
  if (matrix == NULL)
	exit(EXIT_FAILURE);
  return matrix;
}

char* alloc_string (int dimension)
{
  char* string = (char *)malloc((dimension+1)*sizeof(char));
  if (string == NULL)
    exit(EXIT_FAILURE);
  string[dimension]='\0';
  return string;
}

char* alloc_tab (int dim1,int dim2)
{
  int i;
  char* string = (char *)malloc((dim1)*(dim2+1)*sizeof(char));
  if (string == NULL)
    exit(EXIT_FAILURE);
  for (i=0;i<dim1;i++)
    {
      string[ij2k(i,dim2-1,dim2)]='\0';
     }
   return string;
}

/*free memory*/
void free_vector(int * v)
{
  free(v);
}

void free_vector_d(double * v)
{
  free(v);
}

void free_string(char * string)  

{
  free(string);
}


void free_node(Node node)
{
  free_vector_d(node->sol_x);
  
}
/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL */
/*A passer dans utilities */
void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) 
    {
      free (*ptr);
      *ptr = NULL;
    }
}  

/*utilities*/
void copy_line(double * B, double * v, int i,int t)
{
  int j;
  for(j=0;j<t;j++)
    B[ij2k(i,j,t)] = v[j];
} 

int * copy_vector(int *v, int n)
   {
     int i;
     int * c = alloc_vector(n);
     for (i=0;i<n;i++)
       c[i] = v[i];
     return c;
   }

double * copy_vector_d(double *v, int n)
   {
     int i;
     double * c = alloc_vector_d(n);
     if(v!=NULL)
       {
	 for (i=0;i<n;i++)
	   c[i] = v[i];
       }
     else
       {
       for (i=0;i<n;i++)
	 c[i] = 0;
       }
     return c;
   }


double * copy_matrix_d(double *matrix, int m,int n)
   {
     int i;
     
     double* new_matrix = (double *)malloc(m*n*sizeof(double));
     for (i=0;i<m*n;i++)
       new_matrix[i] = matrix[i];
     return new_matrix;
   }

double * copy_matrix_3d(double *matrix, int m,int n, int p)
   {
     int i;
     
     double* new_matrix = (double *)malloc(m*n*p*sizeof(double));
     for (i=0;i<m*n*p;i++)
       new_matrix[i] = matrix[i];
     return new_matrix;
   }

char * copy_string(char *s)
   {
     int i =0;
     char * a = alloc_string(25);
     while (s[i] != '\0')
       {
	 a[i] = s[i];
	 i++;
       }
     a[i]='\0';
     return a;
   }

double v_abs(double a)
{
  if (Positive_BB(a))
    return a;
  return -a;

}

void v_abs_ref(double *a)
 {
  if (Negative_BB(*a))
    *a = -(*a);
} 

int sum_vector(int * a, int j,int n)
{
  int i;
  int b=0;
  for(i=j;i<j+n;i++)
    b+=a[i];
  return b;
}

double sum_vector_d(double * a, int j,int n)
{
  int i;
  double b=0;
  for(i=j;i<j+n;i++)
    b+=a[i];
  return b;
}

double * square_vector(double * v, int n)
{
  int i;
  double * r =alloc_vector_d(n);
   for(i=0;i<n;i++)
     r[i]=v[i]*v[i];
   return r;
}

int random_ind (int a, int b)
{

  return (rand() % (b-a) +a);
} 

double random_double(double min, double max) 
{
    //used it to generate new number each time
    srand( (unsigned int) time(NULL) );

    double randomNumber, range , tempRan, finalRan;
    //generate random number form 0.0 to 1.0
    randomNumber = (double)rand() / (double)RAND_MAX;
    //total range number from min to max eg. from -2 to 2 is 4
    //range used it to pivot form -2 to 2 -> 0 to 4 for next step
    range = max - min;
    //illustrate randomNumber to range
    //lets say that rand() generate 0.5 number, thats it the half 
    //of 0.0 to 1.0, show multiple range with randomNumber we get the
    //half in range. eg 4 * 0.5 = 2
    tempRan = randomNumber * range;
    //add the min to tempRan to get the correct random in ours range
    //so in ours example we have: 2 + (-2) = 0, thats the half in -2 to 2
    finalRan = tempRan + min;
    return finalRan;
}


int maximum(int * vector,int n)
{
  int i;
  int max;
  max = vector[0];  
  for(i=1;i<n;i++)
    if(vector[i] > max)
      max = vector[i];
  return max;
}

double maximum_d(double * vector,int n)
{
  int i;
  double max;
  max = vector[0];  
  for(i=1;i<n;i++)
    if(vector[i] > max)
	max = vector[i];
  
  return max;
}

void minimum_d_0(double * vector,int n, int * ind, double *min)
{
  int i;
  *ind=0;
  *min = vector[0];  
  for(i=1;i<n;i++)
    if(*min > vector[i])
      {
	*min = vector[i];
	*ind=i;
      }
  vector[*ind]=0;
}

void maximum_d_0(double * vector,int n, int * ind, double *max)
{
  int i;
  *ind=0;
  *max = vector[0];  
  for(i=1;i<n;i++)
      if(vector[i]>*max)
	{
	  *max = vector[i];
	  *ind=i;
	}
  vector[*ind]=0;
}

double maximum_d_abs(double * vector,int n)
{
  int i;
  double max;
  double * v=(double*)malloc(n*sizeof(double));
  for(i=0;i<n;i++)
    if(Negative_BB(vector[i]))
      v[i] = -vector[i];
    else
      v[i]=vector[i];
  max=maximum_d(v,n);
  return max;
}

double minimum_d(double * vector,int n)
{
  int i;
  double min;
  min = vector[0];
  for(i=1;i<n;i++)
    if(min > vector[i])
	min = vector[i];
  return min;
}

int minimum(int * vector,int n)
{
  int i;
  int min;
  min = vector[0];
  for(i=1;i<n;i++)
      if(min > vector[i])
	min = vector[i];
  return min;
}

double * sum_two_vector(double * v1, int * v2,int n)
{
  int i;
  for(i=1;i<n;i++)
    v1[i]=v1[i]+(double)v2[i];
  return v1;
}


int nb_non_zero_vector(double * v, int n)
{
  int i;
  int c = 0;
  for(i=0;i<n;i++)
    if(!Zero_BB(v[i]))
      c++;
  return c;
} 

int nb_non_zero_vector_int(int * v, int n)
{
  int i;
  int c = 0;
  for(i=0;i<n;i++)
    if(!Zero_BB(v[i]))
      c++;
  return c;
} 

int nb_non_zero_matrix_sym(double * a, int n)
{
  int i,j;
  int c = 0;
  for(i=0;i<n;i++)
    for(j=i;j<n;j++)
      if(!Zero_BB(a[ij2k(i,j,n)]))
	c++;
  return c;
} 

int nb_non_zero_matrix(double *a, int n, int m)
{
  int i,j;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if(!Zero_BB(a[ij2k(i,j,n)]))
	c++;
  return c;
} 

int nb_non_zero_matrix_i(int *a, int n, int m)
{
  int i,j;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if(!Zero_BB(a[ij2k(i,j,n)]))
	c++;
  return c;
} 

int nb_non_zero_matrix_sup_i(int *a, int n, int m)
{
  int i,j;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=i;j<n;j++)
      if(!Zero_BB(a[ij2k(i,j,n)]))
	c++;
  return c;
} 

int nb_non_zero_matrix_without_diag(double *a, int n, int m)
{
  int i,j;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if (i!=j)
	if(!Zero_BB(a[ij2k(i,j,n)]))
	  c++;
  return c;
} 

int nb_non_zero_matrix_without_diag_i(int *a, int n, int m)
{
  int i,j;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      if (i!=j)
	if(!Zero_BB(a[ij2k(i,j,n)]))
	  c++;
  return c;
}

int nb_non_zero_matrix_3d(double * a, int m, int n, int p)
{
  int i,j,k;
  int c = 0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      for(k=0;k<p;k++)
	if(!Zero_BB(a[ijk2l(i,j,k,n,p)]))
	  c++;
  return c;
} 

int is_in_vector(int val,int *v,int n)
{
  int i;
  for(i=0;i<n;i++)
    if (val == v[i])
      return 1;
  return 0;
}

/*computing logarithme in basis 2*/
int log_base_2 (int x)
{
  return (log(x)/log(2));
}


int vector_length_base_2(double * u, int n)
{
  int i,borne_u;
  int N = 0;
  for(i=0;i<n;i++)
    {
      borne_u= (int)log_base_2((int)u[i])+1;
      N = N + borne_u;
    }
  return N;
}

int * create_u_base2(double * u, int n)
{
  int * lgbis =alloc_vector(n);
  int i;
  for(i=0;i<n;i++)
    lgbis[i]= log_base_2((int)u[i]) + 1;
  return lgbis;
}


int * create_u_base2_cumulated(double * u, int n)
{
  int * lg = alloc_vector(n);
  int i;
  lg[0]= log_base_2((int)u[0]) + 1;
  for(i=1;i<n;i++)
    lg[i] = lg[i-1] + log_base_2((int)u[i]) + 1;
  return lg;
}

char * intochar(int val)
{
  int tmp1;
  int i=0;
  int tmp=val;
  int k=0;
  /*i=log base 10 de val*/
  
  while(tmp>=10)
    {
      tmp=(int)(tmp/10);
      i++;
    }
  char *s = alloc_string(i+1);
 
  tmp=val;
  while(i>=0)
    {
      tmp1=(int)(tmp/pow(10,i));
      s[k]=48+tmp1;
      tmp=tmp-tmp1*pow(10,i);
      i--;
      k++;
    }
  
  return s;
}


/*operations on  matrices*/

int is_symetric(double * a, int n)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      if(a[ij2k(i,j,n)] != a[ij2k(j,i,n)])
	return 0;
  return 1;
}
int is_symetric_3d(double * a,int k, int n)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      if(a[ijk2l(k,i,j,n,n)] != a[ijk2l(k,j,i,n,n)])
	return 0;
  return 1;
}
double * transpose_matrix( double * a, int n)
{
  int i,j;
  double * r = alloc_matrix_d(n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      r[ij2k(i,j,n)] = a[ij2k(j,i,n)];

  return r;
} 

double * transpose_matrix_mn( double * a, int m, int n)
{
  int i,j;
  double * r = alloc_matrix_d(n,m);
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      r[ij2k(j,i,m)] = a[ij2k(i,j,n)];

  return r;
} 

double * multiplication_matrices( double * a, double * b, int n)
{
  int i,j,k;
  double * r = alloc_matrix_d(n,n);
  double s=0;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	for(k=0;k<n;k++)
	  s= s + a[ij2k(i,k,n)] * b[ij2k(k,j,n)];
      	r[ij2k(i,j,n)] = s;
	s=0;
      }
  return r;
} 

double * opposite_matrix_3d(double * M, int m, int n, int p)
{
  int i,j,k;
  double * new_m = alloc_matrix_3d(m,n,p);
  for (i=0;i<m;i++)
    for(j=0;j<n;j++)
       for(k=0;k<p;k++)
	 new_m[ijk2l(i,j,k,n,p)] = -M[ijk2l(i,j,k,n,p)];
  return new_m;

}
/*******************************************************************************************/
/*********************** utilities for the sdp solver***************************************/
/*******************************************************************************************/

int nb_cont_Xij(int n, int nb_int, int p)
{
  int i,cont = 0;
  for(i=n+p-1;i>=n - nb_int;i--)
    cont=cont+i;
  return cont;
}

int nb_max_cont(SDP psdp, int * ind_X1,int * ind_X2,int *ind_X3,int*ind_X4)
{

  int i,j;
  int cont = nb_cont_Xij(psdp->n,psdp->nb_int,0);

  /* Compute the number of constraints that ensure dual feasibility that are already in the sdp */
 
  
  *ind_X1=1;
  *ind_X2=*ind_X1+cont; 
  *ind_X3=*ind_X2+cont;
  *ind_X4=*ind_X3+cont;
  
  return *ind_X4+ cont-1; 
}

void make_matrix_with_vector(SDP psdp, double ** m)
{
  int i,j;
  int k=0;
  
  m[0]=alloc_matrix_d (psdp->n+1,psdp->n+1);
  if (m[0] == NULL)
    exit(EXIT_FAILURE);
  
  for (i=0;i<1+psdp->n;i++)
    for (j=i;j<1+psdp->n;j++)
      {
	if (i==j)
	  m[0][ij2k(i,j,psdp->n+1)]=psdp->x[k];
	
	if (i !=j)
	  {
	    m[0][ij2k(i,j,psdp->n+1)]=psdp->x[k];
	    m[0][ij2k(j,i,psdp->n+1)]=m[0][ij2k(i,j,psdp->n+1)];
	  }
	k++;
      }
}

void make_matrix_with_vector_new_sub(SDP psdp, double ** m, double * X)
{
  int i,j;
  int k=0;
  
  m[0]=alloc_matrix_d (psdp->n+1,psdp->n+1);
  if (m[0] == NULL)
    exit(EXIT_FAILURE);

  for (i=0;i<1+psdp->n;i++)
    for (j=i;j<1+psdp->n;j++)
      {
	if (i==j)
	  m[0][ij2k(i,j,psdp->n+1)]=X[k];
	
	if (i !=j)
	  {
	    m[0][ij2k(i,j,psdp->n+1)]=X[k];
	    m[0][ij2k(j,i,psdp->n+1)]=m[0][ij2k(i,j,psdp->n+1)];
	  }
	k++;
      }
}

void dualized_objective_function(SDP psdp, double ** q_beta, double ** c_beta, double * l_beta)
{

  int i,j,k;

  *q_beta = copy_matrix_d(psdp->q,psdp->n,psdp->n);  
  *c_beta = copy_vector_d(psdp->c,psdp->n);
  *l_beta=psdp->cons;

  /* for all i,j q_ij = - q_ij and c_i = - c_i */
  for(i=0;i<psdp->n;i++)
    for(j=0;j<psdp->n;j++)
      q_beta[0][ij2k(i,j,psdp->n)] = -q_beta[0][ij2k(i,j,psdp->n)];
  
  for(i=0;i<psdp->n;i++)
    c_beta[0][i] = -c_beta[0][i];

  for(k=0;k<psdp->nb_cont;k++)
    {
      if (psdp->constraints[k] >0)
	{
	  i = psdp->tab[psdp->constraints[k]-1][0];
	  j = psdp->tab[psdp->constraints[k]-1][1];
	  switch (psdp->tab[psdp->constraints[k]-1][2])
	    {
	    case 1:
	      /*Dualization of constraints  x_ix_j - u_ix_j - l_jx_i + u_il_j<= 0*/
	      q_beta[0][ij2k(i,j,psdp->n)]= q_beta[0][ij2k(i,j,psdp->n)] - psdp->beta_1[ij2k(i,j,psdp->n)];
	      q_beta[0][ij2k(j,i,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)];
	      c_beta[0][j] = c_beta[0][j] + 2*psdp->beta_1[ij2k(i,j,psdp->n)]*psdp->u[i];
	      c_beta[0][i] = c_beta[0][i] + 2*psdp->beta_1[ij2k(i,j,psdp->n)]*psdp->l[j];
	      l_beta[0]=l_beta[0] - 2*psdp->beta_1[ij2k(i,j,psdp->n)]*psdp->u[i]*psdp->l[j];
	      break;
	      
	    case 2:
	      /*Dualization of constraints x_ix_j - u_jx_i -l_ix_j + u_jl_i<= 0*/
	      q_beta[0][ij2k(i,j,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)] - psdp->beta_2[ij2k(i,j,psdp->n)];
	      q_beta[0][ij2k(j,i,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)];
	      c_beta[0][i]=c_beta[0][i] + 2*psdp->beta_2[ij2k(i,j,psdp->n)]*psdp->u[j];
	      c_beta[0][j]=c_beta[0][j] + 2*psdp->beta_2[ij2k(i,j,psdp->n)]*psdp->l[i];
	      l_beta[0]=l_beta[0] - 2*psdp->beta_2[ij2k(i,j,psdp->n)]*psdp->u[j]*psdp->l[i];
	      break;
	      
	    case 3:
	      /* Dualization of constraints - x_ix_j + u_jx_i + u_ix_j - u_iu_j <= 0*/
	      q_beta[0][ij2k(i,j,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)] + psdp->beta_3[ij2k(i,j,psdp->n)];
	      q_beta[0][ij2k(j,i,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)];
	      c_beta[0][i]=c_beta[0][i]- 2*psdp->beta_3[ij2k(i,j,psdp->n)]*psdp->u[j];
	      c_beta[0][j]=c_beta[0][j]- 2*psdp->beta_3[ij2k(i,j,psdp->n)]*psdp->u[i];
	      l_beta[0]=l_beta[0] + 2*psdp->beta_3[ij2k(i,j,psdp->n)]*psdp->u[i]*psdp->u[j];
	      break;
	      
	    case 4:
	      /* Dualization of constraints - x_ix_j+ l_jx_i + l_ix_j - l_il_j <= 0*/
	      q_beta[0][ij2k(i,j,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)] + psdp->beta_4[ij2k(i,j,psdp->n)];
	      q_beta[0][ij2k(j,i,psdp->n)]=q_beta[0][ij2k(i,j,psdp->n)];
	      c_beta[0][i]=c_beta[0][i]- 2*psdp->beta_4[ij2k(i,j,psdp->n)]*psdp->l[j];
	      c_beta[0][j]=c_beta[0][j]- 2*psdp->beta_4[ij2k(i,j,psdp->n)]*psdp->l[i];
	      l_beta[0]=l_beta[0] + 2*psdp->beta_4[ij2k(i,j,psdp->n)]*psdp->l[i]*psdp->l[j];
	      break;
	      
	    default:
	      printf("\n Error of constraint type (dualize objective function)\n");
	      break;
	    }
	}
       else
	 break;
    }

  
}


void evaluate_objective(double * objective_value,SDP psdp, double ** q_beta, double ** c_beta, double *l_beta)
{
  /* Objective function value : x^t(- Q - beta^1 - beta^2 + beta^3 + beta^4)x + (- c + 2*beta_1u_i + 2*beta^2u_j - 2*beta^3(u_i+u_j))^tx + 2*beta^3u_iu_j  */

  int i,j;
  
  double * matrix;
  *objective_value=0;
  
  make_matrix_with_vector(psdp, &matrix);
  
  /*f(x) =f(x) + x^t( - Q - beta^1 - beta^2 + beta^3 + beta^4)x*/
  for(i=0;i<psdp->n;i++)
    for(j=0;j<psdp->n;j++)
      *objective_value=*objective_value + matrix[ij2k(i+1,j+1,psdp->n+1)]*q_beta[0][ij2k(i,j,psdp->n)];
  
  /*f(x) =f(x) + (- c + beta_1u_i + beta^2u_j - beta^3(u_i+u_j))^tx*/
  for(i=0;i<psdp->n;i++)
    *objective_value=*objective_value + matrix[ij2k(0,i+1,psdp->n+1)]*c_beta[0][i];
  
  /*f(x) =f(x) + beta^3u_iu_j */
  *objective_value=*objective_value + *l_beta;
 
  free_vector_d(matrix);
}


int reduce(double *m, int n, int a, int b, float factor){
	int i;
	if(m == NULL)
		return 0;
	if(n < a || n < b)
		return 0;
	for(i = 0; i < n; i++){
	  m[ij2k(i,b,n)]  -= m[ij2k(i,a,n)]*factor;
	}

	return 1;
}



void compute_subgradient(double ** check_violated,SDP psdp)
{
  
  double * matrix;
  
  make_matrix_with_vector(psdp, &matrix);
  
  /* Attention les indices de la matrice matrix sont décalés de 1*/

  int i,j,k;
  
 
  for(k=0;k<psdp->nb_cont;k++)
    {
      if (psdp->constraints[k] >0)
	{
	  i = psdp->tab[psdp->constraints[k]-1][0];
	  j = psdp->tab[psdp->constraints[k]-1][1];
	  switch (psdp->tab[psdp->constraints[k]-1][2])
	    {
	    case 1:
	      /* subgradient =  - x_ix_j + u_ix_j + l_jx_i - u_il_j*/ 
	      check_violated[0][k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->u[i]*psdp->l[j];
	      break;
	      
	    case 2:
	      /* subgradient = -x_ix_j + u_jx_i+ l_ix_j - u_jl_i*/
	      check_violated[0][k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)]+ psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)]- psdp->u[j]*psdp->l[i];
	      break;
	      
	    case 3:
	      /* subgradient =  x_ix_j - u_jx_i - u_ix_j + u_iu_j */
	      check_violated[0][k] = matrix[ij2k(i+1,j+1,psdp->n+1)] - psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->u[i]*psdp->u[j];
	      break;
	      
	    case 4:
	      /* subgradient = x_ix_j- l_jx_i - l_ix_j + l_il_j*/
	      check_violated[0][k] = matrix[ij2k(i+1,j+1,psdp->n+1)]- psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->l[i]*psdp->l[j];
	      break;
	      
	    default:
	      printf("\n Erreur de type de contrainte (compute subgradient)\n");
	      break;
	    }
	}
      else
	break;
     
    }

  
  
  free_vector_d(matrix);
}


void compute_new_subgradient(double ** check_violated,SDP psdp,int new_length, int ** variable_indices, double ** X)
{
  
  double * matrix;
  
  make_matrix_with_vector_new_sub(psdp, &matrix,X[0]);
  

/* Attention les indices de la matrice matrix sont décalés de 1*/

  int i,j,k;
  
  int * constraints_new=alloc_vector(new_length);
  int index=variable_indices[0][0];
  
  for(i=0;i<new_length;i++)
    {
      constraints_new[i]=psdp->constraints[index+i];
    }

   for(k=0;k<new_length;k++)
    {
      if (constraints_new[k] >0)
	{
	  i = psdp->tab[constraints_new[k]-1][0];
	  j = psdp->tab[constraints_new[k]-1][1];
	  switch (psdp->tab[constraints_new[k]-1][2])
	    {
	    case 1:
	      /*  new subgradient =  - x_ix_j + u_ix_j  + l_jx_i - u_il_j*/ 
	      check_violated[0][k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)]+ psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->u[i]*psdp->l[j];
	      break;
	      
	    case 2:
	      /* new subgradient = -x_ix_j + u_jx_i+ l_ix_j - u_jl_i*/
	      check_violated[0][k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)]+ psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)]- psdp->u[j]*psdp->l[i];
	      break;
	      
	    case 3:
	      /* new subgradient =  x_ix_j - u_jx_i - u_ix_j + u_iu_j */
	      check_violated[0][k] = matrix[ij2k(i+1,j+1,psdp->n+1)] - psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->u[i]*psdp->u[j];
	      break;
	      
	    case 4:
	      /* new subgradient = x_ix_j- l_jx_i - l_ix_j + l_il_j*/
	      check_violated[0][k] = matrix[ij2k(i+1,j+1,psdp->n+1)]- psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->l[i]*psdp->l[j];
	      break;
	      
	    default:
	      printf("\n Erreur de type de contrainte (compute new subgradient)\n");
	      break;
	      
	    }
	}
      else
	break;
   
    }
  
   free_vector_d(matrix);
}

/* purge or add dualized constraints to the objective function */
int check_constraints_violated(SDP psdp, double ** check_violated, double ** x)
{
  
  int i,j,k;
  double * zcheck_violated=alloc_vector_d(psdp->length);;
  double * matrix;

  free_vector_d(check_violated[0]);
  make_matrix_with_vector_new_sub(psdp, &matrix,x[0]);
  
  for(k=0;k<psdp->length;k++)
    {
      i = psdp->tab[k][0];
      j = psdp->tab[k][1];
     
      switch (psdp->tab[k][2])
	{
	case 1:
	/* check_violated =  - x_ix_j + u_ix_j  + l_jx_i - u_il_j*/ 
	  zcheck_violated[k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)]+ psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)]- psdp->u[i]*psdp->l[j];
	  break;

	case 2:
	   /* check_violated = -x_ix_j + u_jx_i+ l_ix_j- u_jl_i*/
	  zcheck_violated[k] = -matrix[ij2k(i+1,j+1,psdp->n+1)] + psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)]+ psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)]- psdp->u[j]*psdp->l[i];
	  break;
	  
	case 3:
	  /* check_violated =  x_ix_j - u_jx_i - u_ix_j + u_iu_j */
	  zcheck_violated[k] = matrix[ij2k(i+1,j+1,psdp->n+1)] - psdp->u[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->u[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->u[i]*psdp->u[j];
	  break;
	  
	case 4:
	  /* check_violated = x_ix_j- l_jx_i - l_ix_j + l_il_j*/
	  zcheck_violated[k] = matrix[ij2k(i+1,j+1,psdp->n+1)]- psdp->l[j]*matrix[ij2k(0,i+1,psdp->n+1)] - psdp->l[i]*matrix[ij2k(0,j+1,psdp->n+1)] + psdp->l[i]*psdp->l[j];
	  break;

	default:
	  printf("\n Erreur de type de contrainte (check violated) \n");
	  break;

	}

    }

  
  
  check_violated[0]=zcheck_violated;
 
  free_vector_d(matrix);

  for(i=0;i<psdp->length;i++)
    if (Negative_BB(check_violated[0][i]))
      return 0;
  return 1;
}


void purge_constraints_not_violated(SDP psdp, double ** slacks, int * n_assign, int ** assign_new_from_old)
{
  int i;
   
  n_assign[0]=0;
 
  int * zassign_new_from_old=alloc_vector(psdp->nb_cont);
  int * zconstraints=alloc_vector(psdp->nb_cont);
  
  
  /* keep constraints that are still violated*/
  for(i=0;i<psdp->nb_cont;i++)
    {
      if ( Zero_VIOL_CB(slacks[0][i]))
	{
	  zassign_new_from_old[n_assign[0]]=i;
	  zconstraints[n_assign[0]]=psdp->constraints[i];
	  n_assign[0]++;
	}
    }
  for(i=0;i<n_assign[0];i++)
    {
      psdp->constraints[i]=zconstraints[i];
      psdp->constraints_old[i]=zconstraints[i];
    }
  
  for(i;i<psdp->nb_cont;i++)
    {
      psdp->constraints[i]=0;
      psdp->constraints_old[i]=0;
    }
  
  free_vector(*assign_new_from_old);
  psdp->nb_cont=n_assign[0];
  *assign_new_from_old=alloc_vector(*n_assign);
  for(i=0;i<*n_assign;i++)
    {
      assign_new_from_old[0][i]=zassign_new_from_old[i];
    }

  
  free_vector(zassign_new_from_old);
  free_vector(zconstraints);
}



void add_constraints_violated(SDP psdp, double ** check_violated,int * n_append)
{
  int i;
  int ind_min;
  double min;
  n_append[0]=0;


  
  for(i=0;i<psdp->nb_cont;i++)
    if(psdp->constraints[i]!=0)
      check_violated[0][psdp->constraints[i]-1]=0;

 
   
  /*add constraints that are the most violated*/
  for(i=psdp->nb_cont;i<psdp->nb_max_cont;i++)
    {
      minimum_d_0(check_violated[0],psdp->length,&ind_min,&min);
      if(Negative_VIOL_CB(min)){
	psdp->constraints[i]=ind_min+1;
	n_append[0]++;
      }
      else
	break;
    }
  
  

  
  psdp->nb_cont= psdp->nb_cont+ n_append[0];
  
  
}

void initialize_constraints_violated(SDP psdp, double ** check_violated, int * new_length)
{
  
 int i=0;
  int nb_cont_violated = 0;
  TAB_IND * tab =NULL; 
 
  tab = create_tab_ind_from_check_violated(check_violated[0],psdp->length);

 

  quickSort(&tab, 0, psdp->length-1);
    
  
  NB_MAX_CONT = psdp->nb_max_cont;
 
  
  nb_cont_violated = initialize_constraints_tab(&tab,psdp->constraints,psdp->constraints_old,psdp->nb_max_cont);

 
  
  for (i=nb_cont_violated; i<psdp->nb_cont;i++)
    psdp->constraints[i]=0;
 
  
  if (nb_cont_violated < psdp->nb_cont) 
    *new_length=nb_cont_violated;
  else  
    *new_length=psdp->nb_cont;
   
  
}


void update_beta_value(SDP psdp, double * zbeta)
{
  int i,j,k;
    
  /*Attention : on divise par 2 ici car on recopie symetriquement dans la matrice du coup il est necessaire de multiplier par 2 les coefficients des termes lineaires et constants*/
  
  free_vector_d(psdp->beta_1);
  free_vector_d(psdp->beta_2);
  free_vector_d(psdp->beta_3);
  free_vector_d(psdp->beta_4);

  double * zbeta_1=alloc_matrix_d(psdp->n,psdp->n);
  double * zbeta_2=alloc_matrix_d(psdp->n,psdp->n);
  double * zbeta_3=alloc_matrix_d(psdp->n,psdp->n);
  double * zbeta_4=alloc_matrix_d(psdp->n,psdp->n);

  for(i=0;i<psdp->n;i++)
    for(j=0;j<psdp->n;j++)
      {
	zbeta_1[ij2k(i,j,psdp->n)]=0;
	zbeta_2[ij2k(i,j,psdp->n)]=0;
	zbeta_3[ij2k(i,j,psdp->n)]=0;
	zbeta_4[ij2k(i,j,psdp->n)]=0;
      }

  for(k=0;k<psdp->nb_cont;k++)
    {
      if (psdp->constraints[k] >0)
	{
	  i = psdp->tab[psdp->constraints[k]-1][0];
	  j = psdp->tab[psdp->constraints[k]-1][1];
	  switch (psdp->tab[psdp->constraints[k]-1][2])
	    {
	    case 1:
	      /*Dual variables associated to constraints  x_ix_j - u_ix_j  - l_jx_i + u_il_j <= 0*/
	      zbeta_1[ij2k(i,j,psdp->n)]=zbeta[k]/2;
	      zbeta_1[ij2k(j,i,psdp->n)]=zbeta_1[ij2k(i,j,psdp->n)];
	      break;
	      
	    case 2:
	      /*Dual variables associated to constraints x_ix_j - u_jx_i  - l_ix_j + u_jl_i  <= 0*/
	      zbeta_2[ij2k(i,j,psdp->n)]=zbeta[k]/2;
	      zbeta_2[ij2k(j,i,psdp->n)]=zbeta_2[ij2k(i,j,psdp->n)];
	      break;
	      
	    case 3:
	      /*Dual variables associated to constraints - x_ix_j + u_jx_i + u_ix_j - u_iu_j <= 0*/
	      zbeta_3[ij2k(i,j,psdp->n)]=zbeta[k]/2;
	      zbeta_3[ij2k(j,i,psdp->n)]=zbeta_3[ij2k(i,j,psdp->n)];
	      break;
	      
	    case 4:
	      /*Dual variables associated to constraints - x_ix_j+ l_jx_i + l_ix_j - l_il_j <= 0*/
	      zbeta_4[ij2k(i,j,psdp->n)]=zbeta[k]/2;
	      zbeta_4[ij2k(j,i,psdp->n)]=zbeta_4[ij2k(i,j,psdp->n)];
	      break;
	      
	    default:
	      printf("\n Erreur de type de contrainte (update beta) \n");
	      break;
	    }
	}
      else
	break;

    }

 
  psdp->beta_1 =zbeta_1;
  psdp->beta_2 =zbeta_2;
  psdp->beta_3 =zbeta_3;
  psdp->beta_4 =zbeta_4;
 
}


/*******************************************************************************************/
/*********************** utilities for the sdp solver cont ***************************************/
/*******************************************************************************************/



int nb_max_cont_mixed(SDP psdp, int * ind_X1,int * ind_X2,int *ind_X3,int*ind_X4)
{

  int i,j;
  int cont = (psdp->n-1)*psdp->n/2;

  /* Compute the number of constraints that ensure dual feasibility that are already in the sdp */
 
  
  *ind_X1=1;
  *ind_X2=*ind_X1+cont; 
  *ind_X3=*ind_X2+cont;
  *ind_X4=*ind_X3+cont;
  
  return *ind_X4+ cont-1; 
}



/****************************************************************************************/
/****************** utilities for the branch and bound **********************************/
/****************************************************************************************/



double check_value(double a)
{
  double b;
  if (a - (int)a < EPS_INT)
    {
      b= (int )a;
     }
  else
    {
      if ((int)a - a == 0)
	b= (int )a;
      else
	b = (int )a +1;
    }
    
  return b;
}

int is_integer(double a)
{
  if (a - (int)a < EPS_INT)
    return 1;
  else
    if ((int)a+1 - a < EPS_INT)
      return 1;
  return 0;
}

double sum_r_b_r(double * b,int m)
{
  int r;
  double res=0;
  for(r=0;r<m;r++)
    res= res + b[r]*b[r];
  return res;
}

double sum_r_ari_b_r(double * a, double * b,int m, int n,int i)
{
  int r;
  double res=0;
  for(r=0;r<m;r++)
    res= res + a[ij2k(r,i,n)]*b[r];
  return res;
}

double sum_r_ari_arj(double * a, int m,int n, int i, int j)
{
  int r;
  double res=0;
  for(r=0;r<m;r++)
    res= res + a[ij2k(r,i,n)] * a[ij2k(r,j,n)];
  return res;
}


double sum_ij_qi_qj_x_ij(double * q, double * x, int n)
{
  int i,j;
  double res=0;
  for(i=0;i<n;i++)
      for(j=0;j<n;j++)
    res= res + q[ij2k(i,j,n)] * x[i] * x[j] ;
  return res;
}


double sum_i_ci_x_i(double * c, double * x, int n)
{
  int i;
  double res=0;
  for(i=0;i<n;i++)
    res= res + c[i] * x[i];
  return res;
}

double sum_r_alphaq_bq_r(double * b,double* alphaq,int m)
{
  int r;
  double res=0;
  for(r=0;r<m;r++)
    res= res + alphaq[r] * b[r];
  return res;
}

double sum_ij_dqi_dqj_x_ij(double * m, double * x, int n,int r)
{
  int i,j;
  double res=0;
  for(i=0;i<n;i++)
    res= res + 2*m[ijk2l(r,0,i+1,n+1,n+1)] * x[i] ;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      res= res + m[ijk2l(r,i+1,j+1,n+1,n+1)] * x[i] * x[j] ;
  return res;
}


/*Ne marche que pour les coef positifs*/
void refine_bounds(struct miqcp_bab bab){
  int i,j,k;
  int i_0;
  
  double b,ci_0,delta;
  double grad_inf,grad_sup;

 
  
  /* attention que si contraintes quadratiques (inegalites pour le moment)*/
  for (i_0=0;i_0<cqp->n;i_0++)
    {
      for (k=0;k<cqp->pq;k++)
	{
	  /*calcul du nouveau b*/
	  b=cqp->eq[k];
	  for (i=0;i<cqp->n;i++)
	    if (i!=i_0)
	      {
		if (Positive_BB(cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]))
		  b=b-2*cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]*bab.l[i];
		if  (Negative_BB(cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]))
		  b=b-2*cqp->dq[ijk2l(k,0,i+1,cqp->n+1,cqp->n+1)]*bab.u[i];
		
		for (j=0;j<cqp->n;j++)
		  if (j!=i_0)
		    {    
		      if (Positive_BB(cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]))
			b=b-cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]*bab.l[i]*bab.l[j];
		      if (Negative_BB(cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]))
			b=b-cqp->dq[ijk2l(k,i+1,j+1,cqp->n+1,cqp->n+1)]*bab.u[i]*bab.u[j];
		    }
	      }

	  ci_0=2*cqp->dq[ijk2l(k,0,i_0+1,cqp->n+1,cqp->n+1)];
	  for (j=0;j<cqp->n;j++)
	    if (j!= i_0)
	      {
		if (Positive_BB(cqp->dq[ijk2l(k,i_0+1,j+1,cqp->n+1,cqp->n+1)]))
		  ci_0=ci_0 + 2*cqp->dq[ijk2l(k,i_0+1,j+1,cqp->n+1,cqp->n+1)]*bab.l[j];
		if (Negative_BB(cqp->dq[ijk2l(k,i_0+1,j+1,cqp->n+1,cqp->n+1)]))
		  ci_0=ci_0 + 2*cqp->dq[ijk2l(k,i_0+1,j+1,cqp->n+1,cqp->n+1)]*bab.u[j];
	      }
	 
	  if (!Zero_BB(ci_0) && (Zero_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)])))
	    {
	      if  (Positive_BB(b/ci_0) && (b/ci_0 < bab.u[i_0])) 
		{
		  if (i_0 <cqp->nb_int)
		    bab.u[i_0]=(int) (b/ci_0);
		  else
		    bab.u[i_0]=b/ci_0;
		}
	      
	      if (Negative_BB(b/ci_0) && (b/ci_0 > bab.l[i_0]))
		bab.l[i_0]=0;
	    }
	  else
	    if (!Zero_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]))
	      {
		delta = ci_0*ci_0 + 4*b*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)];
		if (Positive_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]) && (!Zero_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)])))
		  {
		    if ((sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]) < bab.u[i_0])
		      	{
			  if (i_0 <cqp->nb_int)
			    bab.u[i_0] = (int)( (sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]));
			  else
			    bab.u[i_0] = (sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]);
			}
		  }
		else
		  if (Negative_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]) && (!Zero_BB(cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)])))
		    {
		      if ((sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]) > bab.l[i_0])
				{
			  if (i_0 <cqp->nb_int)
			    bab.l[i_0] = (int) ((sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]));
			  else
			   bab.l[i_0] = (sqrt(delta) - ci_0)/(2*cqp->dq[ijk2l(k,i_0+1,i_0+1,cqp->n+1,cqp->n+1)]);
			}
		    }
	      }
	  
	  
	
	}
    }

  

  /* pb sans contraintes ou avec ineg quad pour grad_inf*/
  if ((Zero_BB(cqp->m)) && (Zero_BB(cqp->p)) &&(Zero_BB(cqp->mq)) && (Zero_BB(cqp->pq)))
    {
       for (i=0;i<cqp->n;i++)
	 {
	   grad_inf = cqp->c[i];
	   for (j=0;j<cqp->n;j++)
	     if (Positive_BB(cqp->q[ij2k(i,j,cqp->n)]) && (!Zero_BB(cqp->q[ij2k(i,j,cqp->n)])))
	       grad_inf=grad_inf + 2*cqp->q[ij2k(i,j,cqp->n)]*bab.l[j];
	     else
	       grad_inf=grad_inf + 2*cqp->q[ij2k(i,j,cqp->n)]*bab.u[j];
	   
	   if(Positive_BB(grad_inf))
	     bab.u[i]=bab.l[i];
	   
	   
	   
	   if ((Zero_BB(cqp->m)) && (Zero_BB(cqp->p)) &&(Zero_BB(cqp->mq)) && (Zero_BB(cqp->pq)))
	     {
	       grad_sup = cqp->c[i];
	       {
		 for (j=0;j<cqp->n;j++)
		   if (Positive_BB(cqp->q[ij2k(i,j,cqp->n)]) && (!Zero_BB(cqp->q[ij2k(i,j,cqp->n)])))
		     grad_sup=grad_sup + 2*cqp->q[ij2k(i,j,cqp->n)]*bab.u[j];
		   else
		     grad_sup=grad_sup + 2*cqp->q[ij2k(i,j,cqp->n)]*bab.l[j];
		 
		 if(Negative_BB(grad_sup))
		   bab.l[i]=bab.u[i];
		 
	       }
	     }
	 }
       }

  
  

}




/****************************************************************************************/
/****************** utilities for the branch and bound MIQCR ****************************/
/****************************************************************************************/




int feasible_for_constraints(double * x_y){
  int r,k;
  double fx=0;


  for(k=0;k<cqp->nb_int;k++)
    if (!is_integer(x_y[k]))
      return 0;
  
  for(r=0;r<cqp->mq;r++)
    {
      fx= sum_ij_dqi_dqj_x_ij(cqp->aq, x_y, cqp->n,r);
      if (fx != cqp->bq[r])
	return 0;
      fx=0;
    }
  
  for(r=0;r<cqp->pq;r++)
    {
      fx= sum_ij_dqi_dqj_x_ij(cqp->dq, x_y, cqp->n,r);
      if (fx > cqp->eq[r])
	return 0;
      fx=0;
    }

  

    return 1;
}




int  select_i_j_miqcrq(double * x_y, struct miqcp_bab bab , int *i, int *j)
{
  int k,l;
  int feasible=1;
  double max=EPS_BRANCH;
  double slack_sol, slack_select, sum_select,sum_sol;


  /*Priorite sur la diagonale*/
  if (feasible==1)
    {
      for(k=0;k<bab.n;k++)
	if (!Zero_BB(x_y[k] - bab.u[k]) &&  !Zero_BB(x_y[k] - bab.l[k]) && k !=father)
	  {
	    if (!Zero_BB(cqp->nb_var_y[ij2k(k,k,bab.n)])) 
	      {
		slack_sol =  v_abs(x_y[k]*x_y[k] - x_y[cqp->var[(bab.n)+ ij2k(k,k,bab.n)]]);
		slack_select=slack_sol*(cqp->new_q[ij2k(k,k,cqp->n)] - cqp->q[ij2k(k,k,cqp->n)]);
		if ( !Zero_BB(slack_select) && (v_abs(slack_select) > max))
		  {
		    feasible =0;
		    max=v_abs(slack_select);
		    *i=k;
		    *j=k;
		  }
		else
		  if ((v_abs(slack_sol) > EPS_BRANCH) && ( v_abs(slack_sol) > max))
		    {
		      feasible =0;
		      max=v_abs(slack_sol);
		      *i=k;
		      *j=k;
		    }
	      }
	  }
    }
      

  if (feasible==1)
    {
      max=EPS_BRANCH;	
      for(k=0;k<bab.n;k++)
	if (!Zero_BB(x_y[k] - bab.u[k]) &&  !Zero_BB(x_y[k] - bab.l[k]) && k != father)
	  {
	    sum_sol=0;
	    sum_select=0;
	    for(l=k;l<bab.n;l++)
	      {
		if ((!Zero_BB(cqp->nb_var_y[ij2k(k,l,bab.n)]) || !Zero_BB(cqp->nb_var_y[ij2k(l,k,bab.n)])) && !Zero_BB(x_y[l] - bab.u[l]) &&  !Zero_BB(x_y[l] - bab.l[l]))
		  {
		    slack_sol =  v_abs(x_y[k]*x_y[l] - x_y[cqp->var[(bab.n)+ ij2k(k,l,bab.n)]]);
		    slack_select=slack_sol*(cqp->new_q[ij2k(k,l,cqp->n)] - cqp->q[ij2k(k,l,cqp->n)]);
		    sum_sol = sum_sol + slack_sol;
		    sum_select = sum_select + slack_select;
		  }
		if ( !Zero_BB(sum_select) && ( v_abs(sum_select) > max))
  		  {
  		    feasible =0;
  		    max=v_abs(sum_select);
		    *i=k;
  		    *j=l;
  		  }
		else
		  if ((v_abs(sum_sol) > EPS_BRANCH) && ( v_abs(sum_sol) > max))
		    {
		      feasible =0;
		      max=v_abs(sum_sol);
		      *i=k;
		      *j=l;
		    }
	      }
	  }
    }
  
  if (feasible ==0)
    {
     
      father = *i;
      return 1;
    }

  

    for(k=0;k<cqp->nb_int;k++)
      if (!is_integer(x_y[k]))
  	{
  	  *i=k;
  	  *j=k;
	  
  	  return 1;
  	}
  
  
  
  if (feasible)
    return 0;
  
  
  
}



struct miqcp_bab realloc_mbab(struct miqcp_bab bab, int nb_col)
{
  struct miqcp_bab new_bab;
  new_bab.n=bab.n;
  new_bab.p=bab.p;
  new_bab.pq=bab.pq;
  new_bab.nb_int=bab.nb_int;
    
  new_bab.u=alloc_vector_d(nb_col);
  memcpy(new_bab.u,bab.u,nb_col*sizeof(double));
  
  new_bab.l=alloc_vector_d(nb_col);
  memcpy(new_bab.l,bab.l,nb_col*sizeof(double));

  

  
  return new_bab;
}




/*  update variables if i && j < nb_int   */
/*branch 1 (int * int) */
/* x_i <= v_i -1 */
struct miqcp_bab  update_variables_int_branch_1(struct miqcp_bab bab, double * x_y, int* i,int * j)
{
  int k;
  struct miqcp_bab new_bab = realloc_mbab(bab,cqp->nb_col);

  new_bab.i=*i;
  new_bab.j=*j;
  
  /* x_i <= v_i */
  new_bab.u[new_bab.i]=  (int)(x_y[new_bab.i]);
  
  return new_bab; 
}

/*branch 2 (int * int) */
/* x_i >= v_i +1 */
struct miqcp_bab update_variables_int_branch_2(struct miqcp_bab bab, double * x_y, int *i,int *j)
{
  
  int k;
  struct miqcp_bab new_bab= realloc_mbab(bab,cqp->nb_col);

  new_bab.i=*i;
  new_bab.j=*j;
  
  /* x_i >= v_i */
  new_bab.l[new_bab.i]=  (int)(x_y[new_bab.i]) +1;


  return new_bab; 
}


/*  update variables if i < nb_int  && j>= nb_int  */

/*branch 1 (int * cont) */
/* x_i = v_i */
struct miqcp_bab update_variables_mixed_branch_1(struct miqcp_bab bab, double * x_y, int* i,int *j,  C_MIQCP cqp)
{
  int k;
  struct miqcp_bab new_bab = realloc_mbab(bab,cqp->nb_col);
 
  new_bab.i=*i;
  new_bab.j=*j;
  /* x_i = v_i*/
  new_bab.u[new_bab.i]=  check_value(x_y[new_bab.i]);
  new_bab.l[new_bab.i]=  check_value(x_y[new_bab.i]);
   
 
  return new_bab;
}



int test_bound_miqcp(struct miqcp_bab bab,int nb_col )
{
  int i;
  for(i=0;i<bab.n;i++)
    if (bab.u[i] - bab.l[i] <  -EPS_BB) 
      {
	return 0;
      }
  return 1;
  
}


/*  update variables if i && j > nb_int   */


/*branch 1 (cont * cont)  */
/* x_i <= v_i */
struct miqcp_bab  update_variables_cont_branch_1(struct miqcp_bab bab, double * x_y, int* i,int * j)
{
  int k;
  struct miqcp_bab new_bab = realloc_mbab(bab,cqp->nb_col);

  new_bab.i=*i;
  new_bab.j=*j;
  
  /* x_i <= v_i */
  new_bab.u[new_bab.i]=(1- GAMMA) * v_abs((bab.u[new_bab.i] + bab.l[new_bab.i])/2) + GAMMA*x_y[new_bab.i];

  
  return new_bab; 
}

/*branch 2 (cont * cont) */
/* x_i >= v_i */
struct miqcp_bab update_variables_cont_branch_2(struct miqcp_bab bab, double * x_y, int *i,int *j)
{
  
  int k;
  struct miqcp_bab new_bab= realloc_mbab(bab,cqp->nb_col);

  new_bab.i=*i;
  new_bab.j=*j;
  
  /* x_i >= v_i */
  new_bab.l[new_bab.i] =(1- GAMMA) * v_abs((bab.u[new_bab.i] + bab.l[new_bab.i])/2)+ GAMMA* x_y[new_bab.i];
  
  
  return new_bab; 
}


/*****************************************************************************/
/*************** Utilities for changing bounds ************************/ 
/*****************************************************************************/


/*come back to the original bounds*/
void change_bounds_back(){
  double * zx = copy_vector_d(best_x,qp->n);
  
  int i;
  for(i=0;i<qp->n;i++)
    best_x[i] = zx[i]*(qp->init_u[i] - qp->init_l[i]) + qp->init_l[i];
}


/* Change bounds to 0 and 1*/
void change_bounds(){

  int i,j,k;
  double * zq = copy_matrix_d(qp->q, qp->n,qp->n);
  double * zc = copy_vector_d(qp->c,qp->n);
  double *zb,*ze,*zbq,*zeq,* za,*zd,*zaq,*zdq;

  double * zcq = alloc_vector_d(qp->n);
  double consq=0;
  double cons=0;
  
  
  
  if(!Zero_BB(qp->m))
    {
      za = copy_matrix_d(qp->a, qp->m,qp->n);
      zb = copy_vector_d(qp->b,qp->m);
    }  
  else
    {
      za=alloc_matrix_d(1,1);
      zb = copy_vector_d(qp->b,qp->m);
    }
  
  if(!Zero_BB(qp->p))
    {
      zd= copy_matrix_d(qp->d, qp->p,qp->n);
      ze = copy_vector_d(qp->e,qp->p);
    }
  else
    {
      zd=alloc_matrix_d(1,1);
      ze = copy_vector_d(qp->e,qp->p);
    }
  
  if(!Zero_BB(qp->mq))
    {
      zaq = copy_matrix_3d(qp->aq, qp->mq,qp->n+1,qp->n+1);
      zbq = copy_vector_d(qp->bq,qp->mq);
    }  
  else
    {
      zaq=alloc_matrix_d(1,1);
      zbq = copy_vector_d(qp->bq,qp->mq);
    }
  
  if(!Zero_BB(qp->pq))
    {
      zdq= copy_matrix_3d(qp->dq, qp->pq,qp->n+1,qp->n+1);
      zeq = copy_vector_d(qp->eq,qp->pq);
    }
  else
    {
      zdq=alloc_matrix_d(1,1);
      zeq = copy_vector_d(qp->eq,qp->pq);
    }
  

  
  qp->cons= 0;
  /*Objective function*/
  for(i=0;i<qp->n;i++)
    {
      qp->cons= qp->cons + qp->c[i]*qp->l[i];
      qp->c[i]=qp->c[i]*(qp->u[i] - qp->l[i]);
      for(j=0;j<qp->n;j++)
	{
	  qp->cons= qp->cons + zq[ij2k(i,j,qp->n)]*qp->l[i]*qp->l[j];
	  qp->c[i]=qp->c[i] + 2*zq[ij2k(i,j,qp->n)]*(qp->u[i] - qp->l[i])*qp->l[j]  ;
	  qp->q[ij2k(i,j,qp->n)] = zq[ij2k(i,j,qp->n)]*(qp->u[i] - qp->l[i])*(qp->u[j] - qp->l[j]); 
	}
    }
 
  for(i=0;i<qp->n;i++)
    for(j=i+1;j<qp->n;j++)
      {
	qp->q[ij2k(i,j,qp->n)] = (qp->q[ij2k(i,j,qp->n)] + qp->q[ij2k(i,j,qp->n)])/2 ;
	qp->q[ij2k(j,i,qp->n)] = qp->q[ij2k(i,j,qp->n)] ;
	
      }

  
  /******************  linear constraints   ********************************************/
  /*equality Constraints */
  for(k=0;k<qp->m;k++)
    {
      cons = 0;
      for(i=0;i<qp->n;i++)
	{
	  cons =  cons +  za[ij2k(k,i,qp->n)]*qp->l[i];
	  qp->a[ij2k(k,i,qp->n)] =  za[ij2k(k,i,qp->n)]*( qp->u[i] - qp->l[i]);
	}
      
      qp->b[k] = qp->b[k] - cons;
      
    }

  /*inequality Constraints */
  for(k=0;k<qp->p;k++)
    {
      cons = 0;
      for(i=0;i<qp->n;i++)
	{
	  cons =  cons +  zd[ij2k(k,i,qp->n)]*qp->l[i];
	  qp->d[ij2k(k,i,qp->n)] =  zd[ij2k(k,i,qp->n)]*( qp->u[i] - qp->l[i]);
	}
      
      qp->e[k] = qp->e[k] - cons;
      
    }

  
  /******************  quadratic constraints   ********************************************/
  /*equality Constraints */

   
  for(k=0;k<qp->mq;k++)
    {
      consq = 0;
      for(i=0;i<qp->n;i++)
	zcq[i]=0;

      
      for(i=0;i<qp->n;i++)
	{
	  consq =  consq +  zaq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*qp->l[i]*2;
	  zcq[i]= zcq[i] + zaq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*( qp->u[i] - qp->l[i]);
	  for(j=0;j<qp->n;j++)
	    {
	      consq= consq + zaq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[i]*qp->l[j];
	      zcq[i]=zcq[i] + zaq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*( qp->u[i] - qp->l[i])*qp->l[j] ;
	      qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=zaq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*(qp->u[i] - qp->l[i])*(qp->u[j] - qp->l[j]); 
	    }
	}
      
      qp->bq[k] = qp->bq[k] -consq;
      for(i=0;i<qp->n;i++)
	{
	  qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]= zcq[i];
	  qp->aq[ijk2l(k,i+1,0,qp->n+1,qp->n+1)]=zcq[i];
	}
    }

  for(k=0;k<qp->mq;k++)
    for(i=0;i<qp->n;i++)
      for(j=i+1;j<qp->n;j++)
	{
	  qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=(qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] +  qp->aq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)])/2;
	  qp->aq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)] = qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)];
	}


  /*inequality Constraints */
  for(k=0;k<qp->pq;k++)
    {
      consq=0;
      for(i=0;i<qp->n;i++)
	zcq[i]= 0;
	
	
      
      for(i=0;i<qp->n;i++)
	{
	  consq =  consq +  zdq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*qp->l[i]*2;
	  zcq[i]= zcq[i] + zdq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*( qp->u[i] - qp->l[i]);
	  for(j=0;j<qp->n;j++)
	    {
	      consq= consq + zdq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[i]*qp->l[j];
	      zcq[i]=zcq[i] + zdq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]* (qp->u[i] - qp->l[i])*qp->l[j] ;
	      qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=zdq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*(qp->u[i] - qp->l[i])*(qp->u[j] - qp->l[j]); 
	    }
	}
      
      qp->eq[k] = qp->eq[k] -consq;
      for(i=0;i<qp->n;i++)
	{
	  qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]= zcq[i];
	  qp->dq[ijk2l(k,i+1,0,qp->n+1,qp->n+1)]= zcq[i];
	}
    }

  for(k=0;k<qp->pq;k++)
    for(i=0;i<qp->n;i++)
      {
	for(j=i+1;j<qp->n;j++)
	  {
	    qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=(qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)] +  qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)])/2;
	    qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)] = qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)];
	  }
      }
    
  /*upper and lower bounds*/
  for(i=0;i<qp->n;i++){
    if (qp->u[i] == qp->l[i])
      {
    	qp->u[i]=qp->l[i];
    	qp->l[i]=qp->l[i];
      }
    else
      {
	qp->u[i]=1;
	qp->l[i]=0;
      }
  }
  
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
      {
	if ((k==j) && (qp->u[k] == qp->l[k]))
	  {
	    qp->u[i]=qp->l[k]*qp->l[k];
	    qp->l[i]=qp->l[k]* qp->l[k];
	    i++;
	  }
	else
	  {
	    qp->u[i]=qp->u[j]*qp->u[k];
	    qp->l[i]=qp->l[j]*qp->l[k];
	    i++;
	  }
      }
  free_vector_d(zq);
  free_vector_d(zc);
  free_vector_d(zb);
  free_vector_d(ze);
  free_vector_d(zbq);
  free_vector_d(zeq);
  free_vector_d(za);
  free_vector_d(zd);
  free_vector_d(zaq);
  free_vector_d(zdq);
  free_vector_d(zcq);

}

void change_bounds_gen_int(){

  int i,j,k;
 
  double consq=0;
  double cons=0;
  
  
   double * zcq = alloc_vector_d(qp->n);
 

  qp->cons= 0;
  /*Objective function*/
  for(i=0;i<qp->n;i++)
    {
      qp->cons= qp->cons + qp->c[i]*qp->l[i];
      for(j=0;j<qp->n;j++)
	{
	  qp->c[i] =  qp->c[i] + 2*qp->q[ij2k(i,j,qp->n)]*qp->l[j];
	  qp->cons= qp->cons + qp->q[ij2k(i,j,qp->n)]*qp->l[i]*qp->l[j];
	}
    }

  
  /******************  linear constraints   ********************************************/
  /*equality Constraints */
  for(k=0;k<qp->m;k++)
    {
      cons = 0;
     
      for(i=0;i<qp->n;i++)
	{
	  cons =  cons +  qp->a[ij2k(k,i,qp->n)]*qp->l[i];
	  
	}
      
      qp->b[k] = qp->b[k] - cons;
      
    }
  
  /*inequality Constraints */
  for(k=0;k<qp->p;k++)
    {
      cons = 0;
      for(i=0;i<qp->n;i++)
	{
	  cons =  cons +  qp->d[ij2k(k,i,qp->n)]*qp->l[i];
	}
      
      qp->e[k] = qp->e[k] - cons;
      
      
    }
  
  /******************  quadratic constraints   ********************************************/
  /*equality Constraints */

   
  for(k=0;k<qp->mq;k++)
    {
      consq = 0;
      for(i=0;i<qp->n;i++)
	zcq[i]=0;
      for(i=0;i<qp->n;i++)
	{
	  consq =  consq +  qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*qp->l[i]*2;
	  for(j=0;j<qp->n;j++)
	    {
	      consq= consq + qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[i]*qp->l[j];
	      zcq[i]=zcq[i] + qp->aq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[j] ;
	    }
	}  
      qp->bq[k] = qp->bq[k] -consq;
      for(i=0;i<qp->n;i++)
	{
	  qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]= qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]+zcq[i];
	  qp->aq[ijk2l(k,i+1,0,qp->n+1,qp->n+1)]=qp->aq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]+zcq[i];
	}
    }
    



  /*inequality Constraints */
  for(k=0;k<qp->pq;k++)
    {
      consq=0;
     	
      for(i=0;i<qp->n;i++)
	zcq[i]= 0;

   
      for(i=0;i<qp->n;i++)
	{
	  consq =  consq +  qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]*qp->l[i]*2;
	  for(j=0;j<qp->n;j++)
	    {
	      consq= consq + qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[i]*qp->l[j];
	      zcq[i]=zcq[i] + qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]*qp->l[j] ;
	    }
	}
      
      qp->eq[k] = qp->eq[k] -consq;
      for(i=0;i<qp->n;i++)
	{
	  qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]=qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]+ zcq[i];
	  qp->dq[ijk2l(k,i+1,0,qp->n+1,qp->n+1)]=qp->dq[ijk2l(k,0,i+1,qp->n+1,qp->n+1)]+ zcq[i];
	}
    }

 
 
  /*upper and lower bounds*/
  for(i=0;i<qp->n;i++){
    if (qp->u[i] == qp->l[i])
      {
    	qp->u[i]=qp->l[i];
    	qp->l[i]=qp->l[i];
      }
    else
      {
	qp->u[i]=qp->u[i]-qp->l[i];
	qp->l[i]=0;
      }
  }
  
  for(j=0;j<qp->n;j++)
    for(k=j;k<qp->n;k++)
      {
	if ((k==j) && (qp->u[k] == qp->l[k]))
	  {
	    qp->u[i]=qp->l[k]*qp->l[k];
	    qp->l[i]=qp->l[k]* qp->l[k];
	    i++;
	  }
	else
	  {
	    qp->u[i]=qp->u[j]*qp->u[k];
	    qp->l[i]=qp->l[j]*qp->l[k];
	    i++;
	  }
      }
}



void change_bounds_back_gen_int(){
  double * zx = copy_vector_d(best_x,qp->n);
  
  int i;
  for(i=0;i<qp->n;i++)
    best_x[i] = zx[i] + qp->init_l[i];
}



/**********************************************************************************************************/
/***************************** Computation of Hessian and eigen values  **********************************/
/***********************************************************************************************************/
void compute_new_q_boxqp()
{
  double * a_t;
  double * d_t;
  int i,j,k;

  cqp->new_q = alloc_matrix_d(cqp->n,cqp->n);
  

  /*x^T Aq x = bq */
  if (!Zero_BB(cqp->mq))
    {
      cqp->sum_r_aq_alphaq=alloc_matrix_d(cqp->n,cqp->n);
      for(i=0;i<cqp->n;i++)
	for(j=0;j<cqp->n;j++)
	  {
	    cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)]=0;
	    for(k=0;k<cqp->mq;k++)
	     
	      cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)] = cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)] + cqp->aq[ijk2l(k,i+1,j+1, cqp->n+1,cqp->n+1)]*cqp->alphaq[k];
	
	  }

      cqp->sum_r_bq=alloc_vector_d(cqp->n);
      for(i=0;i<cqp->n;i++)
	{
	  cqp->sum_r_bq[i]=0;
	  for(j=0;j<cqp->mq;j++)
	    cqp->sum_r_bq[i] =  cqp->sum_r_bq[i] + cqp->aq[ijk2l(j,0,i+1,cqp->n+1,cqp->n+1)]* cqp->alphaq[j];
	}
       
    }
   
  if (!Zero_BB(cqp->pq))
    {
      cqp->sum_r_dq_alphabisq=alloc_matrix_d(cqp->n,cqp->n);
      for(i=0;i<cqp->n;i++)
	for(j=0;j<cqp->n;j++)
	  {
	    cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)]=0;
	    for(k=0;k<cqp->pq;k++)
	      cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)] = cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)] + cqp->dq[ijk2l(k,i+1,j+1, cqp->n+1,cqp->n+1)]*cqp->alphabisq[k];
	  }

      cqp->sum_r_eq=alloc_vector_d(cqp->n);
      for(i=0;i<cqp->n;i++)
	{
	  cqp->sum_r_eq[i]=0;
	  for(j=0;j<cqp->pq;j++)
	    cqp->sum_r_eq[i] =  cqp->sum_r_eq[i] + cqp->dq[ijk2l(j,0,i+1,cqp->n+1,cqp->n+1)]* cqp->alphabisq[j];
       
	}
    }
   
  for(i=0;i<cqp->n;i++)
    for(j=0;j<cqp->n;j++)
      { 
	cqp->new_q[ij2k(i,j,cqp->n)] = 2*cqp->q[ij2k(i,j,cqp->n)]+ 2*cqp->beta[ij2k(i,j,cqp->n) ];
		
	if (!Zero_BB(cqp->mq))
	  cqp->new_q[ij2k(i,j,cqp->n)] = cqp->new_q[ij2k(i,j,cqp->n)] + 2*cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)];
	
	if (!Zero_BB(cqp->pq))
	  cqp->new_q[ij2k(i,j,cqp->n)]= cqp->new_q[ij2k(i,j,cqp->n)] + 2*cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)];

      }

  
}


void compute_new_q_miqcr()
{
  double * a_t;
  double * d_t;
  int i,j,k;

  qqp->new_q = alloc_matrix_d(qqp->n,qqp->n);
  

  /*x^T Aq x = bq */
  if (!Zero_BB(qqp->mq))
    {
      qqp->sum_r_aq_alphaq=alloc_matrix_d(qqp->n,qqp->n);
      for(i=0;i<qqp->n;i++)
	for(j=0;j<qqp->n;j++)
	  {
	    qqp->sum_r_aq_alphaq[ij2k(i,j,qqp->n)]=0;
	    for(k=0;k<qqp->mq;k++)
	     
		qqp->sum_r_aq_alphaq[ij2k(i,j,qqp->n)] = qqp->sum_r_aq_alphaq[ij2k(i,j,qqp->n)] + qqp->aq[ijk2l(k,i+1,j+1, qqp->n+1,qqp->n+1)]*qqp->alphaq[k];
	
	  }

      qqp->sum_r_bq=alloc_vector_d(qqp->n);
      for(i=0;i<qqp->n;i++)
	{
	  qqp->sum_r_bq[i]=0;
	  for(j=0;j<qqp->mq;j++)
	    qqp->sum_r_bq[i] =  qqp->sum_r_bq[i] + qqp->aq[ijk2l(j,0,i+1,qqp->n+1,qqp->n+1)]* qqp->alphaq[j];
	}
       
    }
   
  if (!Zero_BB(qqp->pq))
    {
      qqp->sum_r_dq_alphabisq=alloc_matrix_d(qqp->n,qqp->n);
      for(i=0;i<qqp->n;i++)
	for(j=0;j<qqp->n;j++)
	  {
	    qqp->sum_r_dq_alphabisq[ij2k(i,j,qqp->n)]=0;
	    for(k=0;k<qqp->pq;k++)
	      qqp->sum_r_dq_alphabisq[ij2k(i,j,qqp->n)] = qqp->sum_r_dq_alphabisq[ij2k(i,j,qqp->n)] + qqp->dq[ijk2l(k,i+1,j+1, qqp->n+1,qqp->n+1)]*qqp->alphabisq[k];
	  }

      qqp->sum_r_eq=alloc_vector_d(qqp->n);
      for(i=0;i<qqp->n;i++)
	{
	  qqp->sum_r_eq[i]=0;
	  for(j=0;j<qqp->pq;j++)
	    qqp->sum_r_eq[i] =  qqp->sum_r_eq[i] + qqp->dq[ijk2l(j,0,i+1,qqp->n+1,qqp->n+1)]* qqp->alphabisq[j];
       
	}
    }
   
 
  for(i=0;i<qqp->n;i++)
    for(j=0;j<qqp->n;j++)
      { 
	qqp->new_q[ij2k(i,j,qqp->n)] = 2*qqp->q[ij2k(i,j,qqp->n)]+ 2*qqp->beta[ij2k(i,j,qqp->n) ];
	
	if (!Zero_BB(qqp->mq))
	  qqp->new_q[ij2k(i,j,qqp->n)] = qqp->new_q[ij2k(i,j,qqp->n)] + 2*qqp->sum_r_aq_alphaq[ij2k(i,j,qqp->n)];
	
	if (!Zero_BB(qqp->pq))
	  qqp->new_q[ij2k(i,j,qqp->n)]= qqp->new_q[ij2k(i,j,qqp->n)] + 2*qqp->sum_r_dq_alphabisq[ij2k(i,j,qqp->n)];

      }
}



void compute_new_lambda_min_miqcr() 
{
  char       jobz, uplo;
  
  int n = qqp->n;
  int info, lwork;
  int lda = qqp->n;

  double * w;
  double* work;

  double * a = copy_matrix_d(qqp->new_q,qqp->n,qqp->n);
  
  jobz  =  'V';
  uplo  =  'U';
  
  w = (double*)malloc( n*sizeof(double) );
  lwork = 5*n;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Solve eigenproblem */
  dsyev_(&jobz,&uplo, &n, a, &lda, w, work, &lwork,&info );

  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }

  qqp->new_lambda_min = w[0]; 
  if (Negative_LM(qqp->new_lambda_min))
    qqp->new_lambda_min=  qqp->new_lambda_min - EPS_LM;
  else
    qqp->new_lambda_min= - EPS_LM;
}


  


void compute_new_q_sbb_miqcr()
{
  double * a_t;
  double * d_t;
  int i,j,k;

  cqp->new_q = alloc_matrix_d(cqp->n,cqp->n);
  

  /*x^T Aq x = bq */
  if (!Zero_BB(cqp->mq))
    {
      cqp->sum_r_aq_alphaq=alloc_matrix_d(cqp->n,cqp->n);
      for(i=0;i<cqp->n;i++)
	for(j=0;j<cqp->n;j++)
	  {
	    cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)]=0;
	    for(k=0;k<cqp->mq;k++)
	     
		cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)] = cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)] + cqp->aq[ijk2l(k,i+1,j+1, cqp->n+1,cqp->n+1)]*cqp->alphaq[k];
	
	  }

      cqp->sum_r_bq=alloc_vector_d(cqp->n);
      for(i=0;i<cqp->n;i++)
	{
	  cqp->sum_r_bq[i]=0;
	  for(j=0;j<cqp->mq;j++)
	    cqp->sum_r_bq[i] =  cqp->sum_r_bq[i] + cqp->aq[ijk2l(j,0,i+1,cqp->n+1,cqp->n+1)]* cqp->alphaq[j];
	}
       
    }
   
  if (!Zero_BB(cqp->pq))
    {
      cqp->sum_r_dq_alphabisq=alloc_matrix_d(cqp->n,cqp->n);
      for(i=0;i<cqp->n;i++)
	for(j=0;j<cqp->n;j++)
	  {
	    cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)]=0;
	    for(k=0;k<cqp->pq;k++)
	      cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)] = cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)] + cqp->dq[ijk2l(k,i+1,j+1, cqp->n+1,cqp->n+1)]*cqp->alphabisq[k];
	  }

      cqp->sum_r_eq=alloc_vector_d(cqp->n);
      for(i=0;i<cqp->n;i++)
	{
	  cqp->sum_r_eq[i]=0;
	  for(j=0;j<cqp->pq;j++)
	    cqp->sum_r_eq[i] =  cqp->sum_r_eq[i] + cqp->dq[ijk2l(j,0,i+1,cqp->n+1,cqp->n+1)]* cqp->alphabisq[j];
       
	}
    }
  
  for(i=0;i<cqp->n;i++)
    for(j=0;j<cqp->n;j++)
      { 
	cqp->new_q[ij2k(i,j,cqp->n)] = 2*cqp->q[ij2k(i,j,cqp->n)]+ 2*cqp->beta[ij2k(i,j,cqp->n) ];
	
	if (!Zero_BB(cqp->mq))
	  cqp->new_q[ij2k(i,j,cqp->n)] = cqp->new_q[ij2k(i,j,cqp->n)] + 2*cqp->sum_r_aq_alphaq[ij2k(i,j,cqp->n)];
	
	if (!Zero_BB(cqp->pq))
	  cqp->new_q[ij2k(i,j,cqp->n)]= cqp->new_q[ij2k(i,j,cqp->n)] + 2*cqp->sum_r_dq_alphabisq[ij2k(i,j,cqp->n)];

      }

}



void compute_new_lambda_min_sbb_miqcr() 
{
  
  char       jobz, uplo;
  
  int n = cqp->n;
  int info, lwork;
  int lda = cqp->n;

  double * w;
  double* work;

  double * a = copy_matrix_d(cqp->new_q,cqp->n,cqp->n);
  
  jobz  =  'V';
  uplo  =  'U';
  
  w = (double*)malloc( n*sizeof(double) );
  lwork = 5*n;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Solve eigenproblem */
  dsyev_(&jobz,&uplo, &n, a, &lda, w, work, &lwork,&info );

  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }

  cqp->new_lambda_min = w[0]; 
  if (Negative_LM(cqp->new_lambda_min))
    cqp->new_lambda_min=  cqp->new_lambda_min - EPS_LM;
  else
    cqp->new_lambda_min= - EPS_LM;
}

