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

#include "in_out.h"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include"utilities.h"
#include"quad_prog.h"

extern char * buf_in;
extern MIQCP qp;
extern double EPS_BETA;
    

/*********************************************************************/
/**********************  Read the instances **************************/
/*********************************************************************/

double * scan_vector_u(FILE * temp, int t)
{
  int i,ret;
  double * v= alloc_vector_d(t + (t+1)*t/2);
  for(i=0;i<t;i++)
    ret=fscanf(temp,"%lf",&v[i]);
   
  return v;
}

double * scan_vector_d(FILE * temp, int t)
{
  int i,nb_non_zero,k,ret;
  double * v= alloc_vector_d(t);
  for(i=0;i<t;i++)
    v[i]=0;
  ret =fscanf(temp,"%d",&nb_non_zero); 
  for(i=0;i<nb_non_zero;i++) 
    {
      ret = fscanf(temp,"%d",&k);
      ret = fscanf(temp,"%lf",&v[k]);
    }
  return v;
}

double * scan_matrix_d(FILE * temp,int t1, int t2)
{
  int i,j,k,ret;
  long nb_non_zero;
  double * M= alloc_matrix_d(t1,t2);
  for(i=0;i<t1;i++)
    for(j=0;j<t2;j++)
      M[ij2k(i,j,t2)]=0;
  ret = fscanf(temp,"%ld",&nb_non_zero);
 
  for(i=0;i<nb_non_zero;i++) 
    {
      ret = fscanf(temp,"%d",&j);
      ret = fscanf(temp,"%d",&k);
      ret = fscanf(temp,"%lf",&M[ij2k(j,k,t2)]);
    }  
  return M;
}

double * scan_matrix_3d(FILE * temp,int t1, int t2)
{
  int i,j,k,l,ret;
  long nb_non_zero;
  double * M= alloc_matrix_3d(t1,t2,t2);
  for(i=0;i<t1;i++)
    for(j=0;j<t2;j++)
      for(l=0;l<t2;l++)
	M[ijk2l(i,j,l,t2,t2)]=0;
  
  ret = fscanf(temp,"%ld",&nb_non_zero);
 
  for(i=0;i<nb_non_zero;i++) 
    {
      ret =  fscanf(temp,"%d",&j);
      ret = fscanf(temp,"%d",&k);
      ret =  fscanf(temp,"%d",&l);
      ret =  fscanf(temp,"%lf",&M[ijk2l(j,k,l,t2,t2)]);
    }  
  return M;
}

/*Read a standard file*/
int read_file()
{
  int i,j,k,l,ret;
  FILE *file_inst;
  char tmp;
  
  /* Gestion des cont. lin.*/
  double * a_t_a;
  double * a_t;
  double b_q;
  

  double * zaq;
  double * zbq;
  double * zdq;
  double * zeq;
  
  qp =new_miqcp();
  
  file_inst = fopen(buf_in,"r");
  if (file_inst != NULL){
    ret =  fscanf(file_inst,"%d",&qp->n);
    ret =  fscanf(file_inst,"%d",&qp->nb_int);
    ret =  fscanf(file_inst,"%d",&qp->m);
    ret = fscanf(file_inst,"%d",&qp->p);
    ret = fscanf(file_inst,"%d",&qp->mq);
    ret = fscanf(file_inst,"%d",&qp->pq);

    int nbcols=qp->n+(qp->n+1)*qp->n/2;
    /*fscanf(file_inst,"%c",&tmp);*/
    ret =  fscanf(file_inst,"%c",&tmp);
    while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
      ret =  fscanf(file_inst,"%c",&tmp);
    ret = fscanf(file_inst,"%c",&tmp);
    qp->u=scan_vector_u(file_inst,qp->n);

    ret = fscanf(file_inst,"%c",&tmp);
    while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
      ret = fscanf(file_inst,"%c",&tmp);
    ret = fscanf(file_inst,"%c",&tmp);
    qp->l=scan_vector_u(file_inst,qp->n);

    qp->init_l = copy_vector_d(qp->l,qp->n);
    qp->init_u = copy_vector_d(qp->u,qp->n);

    i=qp->n;
    for(j=0;j<qp->n;j++)
      for(k=j;k<qp->n;k++)
     	{
	  if (k==j)
	    {
	      if ((Positive_BB(qp->u[j])||Zero_BB(qp->u[j])) && (Positive_BB(qp->l[j]) || Zero_BB(qp->l[j])))
		{
		  qp->u[i]=qp->u[j]*qp->u[j];
		  qp->l[i]=qp->l[j]*qp->l[j];
		  i++;
		}
	      else
		{
		  if ((Positive_BB(qp->u[j])|| Zero_BB(qp->u[j])) && Negative_BB(qp->l[j]))
		    {
		      if (qp->u[j]*qp->u[j] >qp->l[j]*qp->l[j]) 
			qp->u[i]=qp->u[j]*qp->u[j];
		      else
			qp->u[i]=qp->l[j]*qp->l[j];
		      qp->l[i]=qp->l[j]*qp->u[j];
		      i++;
		    }
		  else
		    if (Negative_BB(qp->u[j]) && Negative_BB(qp->l[j]))
		      {
			qp->u[i]=qp->l[j]*qp->l[j];
			qp->l[i]=qp->u[j]*qp->u[j];
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
			qp->u[i]=qp->u[j]*qp->u[k];
		      else
			qp->u[i]=qp->l[j]*qp->l[k];
		      if  (qp->l[j]*qp->u[k] > qp->l[k]*qp->u[j])
			qp->l[i]=qp->l[j]*qp->u[k];
		      else
			qp->l[i]=qp->l[k]*qp->u[j];
		      i++;
		    }
		  
		  else
		    {
		      if (Negative_BB(qp->l[j]))
			{
			  qp->u[i]=qp->u[j]*qp->u[k];
			  qp->l[i]=qp->l[j]*qp->u[k];
			  i++;
			}
		      else
			{
			  if (Negative_BB(qp->l[k]))
			    {
			      qp->u[i]=qp->u[j]*qp->u[k];
			      qp->l[i]=qp->u[j]*qp->l[k];
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
		}
	      else
		{
		  if (Negative_BB(qp->u[j]) && Negative_BB(qp->u[k]) )
		    {
		      qp->u[i]=qp->l[k]*qp->l[j];
		      qp->l[i]=qp->u[k]*qp->u[j];
		      i++;
		    }
		  else
		    {
		      if (Negative_BB(qp->u[j]))
			{
			  qp->u[i]=qp->l[k]*qp->l[j];
			  qp->l[i]=qp->u[k]*qp->l[j];
			  i++;
			}
		      else
			if (Negative_BB(qp->u[k]))
			  {
			    qp->u[i]=qp->l[k]*qp->l[j];
			    qp->l[i]=qp->u[j]*qp->l[k];
			    i++;
			  }
		      
		    }
		}
	      
	    }

	}


    ret = fscanf(file_inst,"%c",&tmp);
    while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
      ret = fscanf(file_inst,"%c",&tmp);

    ret = fscanf(file_inst,"%c",&tmp);
  
    qp->q = scan_matrix_d(file_inst,qp->n,qp->n);
    /*check if Q is symetric : if not make it symetric */

  

    if (!is_symetric(qp->q,qp->n))
      {
	for(i=0;i<qp->n;i++)
	  {
	    for(j=i+1;j<qp->n;j++)
	      {
		qp->q[ij2k(i,j,qp->n)] = (qp->q[ij2k(i,j,qp->n)] +qp->q[ij2k(j,i,qp->n)])/2; 
		qp->q[ij2k(j,i,qp->n)] = qp->q[ij2k(i,j,qp->n)]; 
	      }
	  }
      }

  
    ret = fscanf(file_inst,"%c",&tmp);
    while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
      ret = fscanf(file_inst,"%c",&tmp);

    ret = fscanf(file_inst,"%c",&tmp);

    qp->c = scan_vector_d(file_inst,qp->n);
  
    qp->cons=0;
    if (qp->flag_const)
      {
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
      
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%lf",&qp->cons);
      }

 
    
    /*Ax=b*/
    if(Positive_BB(qp->m))
      {
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	qp->a= scan_matrix_d(file_inst,qp->m,qp->n);

	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	qp->b= scan_vector_d(file_inst,qp->m);
 
	/* creation of Matrix A^t A*/
	a_t_a = alloc_matrix_d(qp->n,qp->n);
	a_t = transpose_matrix_mn(qp->a, qp->m, qp->n);
	 
	for (i=0;i<qp->n;i++)
	  for (j=0;j<qp->n;j++)
	    a_t_a[ij2k(i,j,qp->n)]=0;
      
           
	for (i=0;i<qp->m;i++)
	  for (j=0;j<qp->n;j++)
	    for(k=j;k<qp->n;k++)
	      a_t_a[ij2k(j,k,qp->n)]=a_t_a[ij2k(j,k,qp->n)] + a_t[ij2k(j,i,qp->m)]*qp->a[ij2k(i,k,qp->n)];
      
	for (j=0;j<qp->n;j++)
	  for(k=j+1;k<qp->n;k++)
	    a_t_a[ij2k(k,j,qp->n)]=a_t_a[ij2k(j,k,qp->n)];
      
	b_q=0;
	for(i=0;i<qp->m;i++)
	  b_q = b_q + qp->b[i]*qp->b[i];
    
    
      }
    else
      {
	qp->a=alloc_matrix_d(1,1);
	qp->b=alloc_vector_d(1);
      }

    /*Dx <= e*/
    if(Positive_BB(qp->p))
      {
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	qp->d= scan_matrix_d(file_inst,qp->p,qp->n);

	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	qp->e= scan_vector_d(file_inst,qp->p);
      
      
      }
    else
      {
	qp->d=alloc_matrix_d(1,1);
	qp->e=alloc_vector_d(1);
      }
  
  
  
  
    /* x^T Aq x = bq */
    if(Positive_BB(qp->mq))
      {
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);

	zaq= scan_matrix_3d(file_inst,qp->mq,qp->n+1);
     
	for(k=0;k<qp->mq;k++)
	  if (!is_symetric_3d(zaq,k,qp->n+1))
	    {
	      for(i=0;i<qp->n+1;i++)
		{
		  for(j=i+1;j<qp->n+1;j++)
		    {
		      zaq[ijk2l(k,i,j,qp->n+1,qp->n+1)] = (zaq[ijk2l(k,i,j,qp->n+1,qp->n+1)] + zaq[ijk2l(k,j,i,qp->n+1,qp->n+1)])/2; 
		      zaq[ijk2l(k,j,i,qp->n+1,qp->n+1)] = zaq[ijk2l(k,i,j,qp->n+1,qp->n+1)];
		    }
		}
	    }

	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	zbq= scan_vector_d(file_inst,qp->mq);
     
	if (Positive_BB(qp->m))
	  {
	    qp->mq++;
	  
	    qp->aq=alloc_matrix_3d(qp->mq,qp->n+1,qp->n+1);
	    qp->bq=alloc_vector_d(qp->mq);
	  
	    for (i=0;i<qp->mq-1;i++)
	      for (j=0;j<qp->n+1;j++)
		for(k=0;k<qp->n+1;k++)
		  qp->aq[ijk2l(i,j,k,qp->n+1,qp->n+1)]=zaq[ijk2l(i,j,k,qp->n+1,qp->n+1)];

	  
	    /* cont lin */
	    i=qp->mq-1;
	    qp->aq[ijk2l(i,0,0,qp->n+1,qp->n+1)]=0;
	    for (j=0;j<qp->n;j++)
	      {
		qp->aq[ijk2l(i,0,j+1,qp->n+1,qp->n+1)]=0;
		qp->aq[ijk2l(i,j+1,0,qp->n+1,qp->n+1)]=0;
		for(k=0;k<qp->n;k++)
		  qp->aq[ijk2l(i,j+1,k+1,qp->n+1,qp->n+1)]=a_t_a[ij2k(j,k,qp->n)];
	      }
	    for (i=0;i<qp->mq-1;i++)
	      qp->bq[i] = zbq[i];
	    i=qp->mq-1;
	    qp->bq[i]= b_q;

	 
	  
	  }
	else
	  {
	    qp->aq=alloc_matrix_3d(qp->mq,qp->n+1,qp->n+1);
	    qp->bq=alloc_vector_d(qp->mq);
	  
	    for (i=0;i<qp->mq;i++)
	      for (j=0;j<qp->n+1;j++)
		for(k=0;k<qp->n+1;k++)
		  qp->aq[ijk2l(i,j,k,qp->n+1,qp->n+1)]=zaq[ijk2l(i,j,k,qp->n+1,qp->n+1)];

	    for (i=0;i<qp->mq;i++)
	      qp->bq[i] = zbq[i];
	  }
	free_vector_d(zaq);
	free_vector_d(zbq);
     
      }
    else
      {
	if (Positive_BB(qp->m))
	  {
	    qp->mq++;

	    qp->aq=alloc_matrix_3d(qp->mq,qp->n+1,qp->n+1);
	    qp->bq=alloc_vector_d(qp->mq);
  
	    /* cont lin */
	    i=0;
	    qp->aq[ijk2l(i,0,0,qp->n+1,qp->n+1)]=0;
	    for (j=0;j<qp->n;j++)
	      {
		qp->aq[ijk2l(i,0,j+1,qp->n+1,qp->n+1)]=0;
		qp->aq[ijk2l(i,j+1,0,qp->n+1,qp->n+1)]=0;
		for(k=0;k<qp->n;k++)
		  qp->aq[ijk2l(i,j+1,k+1,qp->n+1,qp->n+1)]=a_t_a[ij2k(j,k,qp->n)];
	      }
	 
	    i=0;
	    qp->bq[i]= b_q;
	   
	  
	  }
	else
	  {
	    qp->aq=alloc_matrix_3d(1,1,1);
	    qp->bq=alloc_vector_d(1);
	  }
      }

    /* x^T Dq x <= eq */
    if(Positive_BB(qp->pq))
      {
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	/* ret = fscanf(file_inst,"%c",&tmp);
	   ret = fscanf(file_inst,"%c",&tmp);*/
	zdq= scan_matrix_3d(file_inst,qp->pq,qp->n+1);
      
	for(k=0;k<qp->pq;k++)
	  if (!is_symetric_3d(zdq ,k,qp->n+1))
	    {
	    
	      for(i=0;i<qp->n+1;i++)
		{
		  for(j=i+1;j<qp->n+1;j++)
		    {
		      zdq[ijk2l(k,i,j,qp->n+1,qp->n+1)] = (zdq[ijk2l(k,i,j,qp->n+1,qp->n+1)] + zdq[ijk2l(k,j,i,qp->n+1,qp->n+1)])/2; 
		      zdq[ijk2l(k,j,i,qp->n+1,qp->n+1)] = zdq[ijk2l(k,i,j,qp->n+1,qp->n+1)];
		    }
		}
	    }

     
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	while((tmp==' ' ) || (tmp =='\r') || (tmp== '\t'))
	  ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	ret = fscanf(file_inst,"%c",&tmp);
	zeq=scan_vector_d(file_inst,qp->pq);
  
	if (Positive_BB(qp->p))
	  {
	    qp->pq=qp->pq+2*qp->n*qp->p;
	  
	    qp->dq=alloc_matrix_3d(qp->pq,qp->n+1,qp->n+1);
	    qp->eq=alloc_vector_d(qp->pq);
	  
	    for (i=0;i<qp->pq-(2*qp->n*qp->p);i++)
	      for (j=0;j<qp->n+1;j++)
		for(k=0;k<qp->n+1;k++)
		  qp->dq[ijk2l(i,j,k,qp->n+1,qp->n+1)]=zdq[ijk2l(i,j,k,qp->n+1,qp->n+1)];
	  
	    /* cont lin */
	    k=qp->pq-2*qp->n*qp->p;
	    for(l=0;l<qp->p;l++)
	      {
		for(i=0;i<qp->n;i++)
		  {
		    k=k+1;
		    qp->dq[ijk2l(k,0,0,qp->n+1,qp->n+1)]=0;
		    for (j=0;j<qp->n;j++)
		      {
		      
			if (i==j)
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)]- qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)]- qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,j+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)];
			  }
			else
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)]/2;
			    qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)]/2;
			  }
		      }


		 
		  }
	      }

	 

	    k=qp->pq-qp->n*qp->p;
	    for(l=0;l<qp->p;l++)
	      {
		for(i=0;i<qp->n;i++)
		  {
		    k=k+1;
		    qp->dq[ijk2l(k,0,0,qp->n+1,qp->n+1)]=0;
		  
		    for (j=0;j<qp->n;j++)
		      {
		      
			if (i==j)
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)]+ qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)]+ qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,j+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)];
			  }
			else
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)]/2;
			    qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)]/2;
			  }
		      }
		  }
	      }
	

	    for (i=0;i<qp->pq-2*qp->n*qp->p;i++)
	      qp->eq[i] = zeq[i];
	    for(i;i<qp->pq-qp->n*qp->p;i++)
	      for(k=0;k<qp->n;k++)
		{
		  qp->eq[i]= -qp->l[k]*qp->e[j];
		  i++;
		}
	    for(j=0;j<qp->p;j++)
	      for(k=0;k<qp->n;k++)
		{
		  qp->eq[i]= qp->u[k]*qp->e[j];
		  i++;
		}
	  }
	else
	  {
	    qp->dq=alloc_matrix_3d(qp->pq,qp->n+1,qp->n+1);
	    qp->eq=alloc_vector_d(qp->pq);
	  
	    for (i=0;i<qp->pq;i++)
	      for (j=0;j<qp->n+1;j++)
		for(k=0;k<qp->n+1;k++)
		  qp->dq[ijk2l(i,j,k,qp->n+1,qp->n+1)]=zdq[ijk2l(i,j,k,qp->n+1,qp->n+1)];
	    for (i=0;i<qp->pq;i++)
	      qp->eq[i] = zeq[i];
	 

	  }

	free_vector_d(zdq);
	free_vector_d(zeq);
    
      }
    else
      {
	if (Positive_BB(qp->p))
	  {
	    qp->pq=qp->pq+2*qp->n*qp->p;
	  
	    qp->dq=alloc_matrix_3d(qp->pq,qp->n+1,qp->n+1);
	    qp->eq=alloc_vector_d(qp->pq);
	  
	    /* cont lin */
	    k=0;
	    for(l=0;l<qp->p;l++)
	      {
		for(i=0;i<qp->n;i++)
		  {
      		
		    qp->dq[ijk2l(k,0,0,qp->n+1,qp->n+1)]=0;
		    for (j=0;j<qp->n;j++)
		      {
		      
			if (i==j)
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)]- qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)]- qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,j+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)];
			  }
			else
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(-qp->l[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)]/2;
			    qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=qp->d[ij2k(l,j,qp->n)]/2;
			  }
		      }


		    k++;
		  }
	      }

	    k=qp->pq-qp->n*qp->p;
	  
	    for(l=0;l<qp->p;l++)
	      {
		for(i=0;i<qp->n;i++)
		  {
		
		    qp->dq[ijk2l(k,0,0,qp->n+1,qp->n+1)]=0;
		  
		    for (j=0;j<qp->n;j++)
		      {
		      
			if (i==j)
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)]+ qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)]+ qp->e[l])/2;
			    qp->dq[ijk2l(k,j+1,j+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)];
			  }
			else
			  {
			    qp->dq[ijk2l(k,0,j+1,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,0,qp->n+1,qp->n+1)]=(qp->u[i]*qp->d[ij2k(l,j,qp->n)])/2;
			    qp->dq[ijk2l(k,j+1,i+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)]/2;
			    qp->dq[ijk2l(k,i+1,j+1,qp->n+1,qp->n+1)]=-qp->d[ij2k(l,j,qp->n)]/2;
			  }
		      }
		    k++;
		  }
	      }
	

	    i=0;
	    for(j=0;j<qp->p;j++)
	      for(k=0;k<qp->n;k++)
		{
		  qp->eq[i]= -qp->l[k]*qp->e[j];
		  i++;
		}
	    for(j=0;j<qp->p;j++)
	      for(k=0;k<qp->n;k++)
		{
		  qp->eq[i]= qp->u[k]*qp->e[j];
		  i++;
		}
	  }
	else
	  {
	    qp->dq=alloc_matrix_3d(1,1,1);
	    qp->eq=alloc_vector_d(1);
	  }
      }

 
  
 
    qp->local_sol = NULL;
  

    fclose(file_inst);
    return 1;
  }
  else
    {
      printf("\ninstance file not found\n");
      return 0;
    }

}
 
