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

#include "liste.h"

#include "quad_prog.h"
#include "utilities.h"

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>

extern double lb_bb;

Liste new_liste()
{
  return  (Liste)malloc(sizeof(struct liste));
}

Node new_node(double value, double * sol, struct miqcp_bab new_bab, int nb_col)
{
   Node node = (Node)malloc(sizeof(struct node));
   if (node == NULL)
    {
      printf ("\nnew_node : Could not create node.\n");
      return NULL;
    }
 
   node->inf_bound_value = value;
   node->sol_x =sol;

   
   node->bab =realloc_mbab_liste(new_bab,nb_col);
   node->previous = NULL;
   node->next =NULL;
   

   return node;
}

Liste init_liste()
{
  Liste liste = new_liste();
  if (liste == NULL)
    {
      printf ("\ninit_liste : Could not create liste.\n");
      return NULL;
    }

  liste->first =NULL;
  liste->size=0;
  return liste;
}


void insert_liste(Liste liste, double value, double * sol, struct miqcp_bab bab, int nb_col)
{
  Node tmp=liste->first;
  int i=0;
  /*init new node*/
  Node node = new_node(value,sol,bab,nb_col);

  
  if (node == NULL)
    {
      printf ("\n insert :Could not create node.\n");
      return;
    }

   if (liste->size==1)
     {
       if (tmp->inf_bound_value < value)
	 {
	   tmp->next=node;
	   node->previous=tmp;
	   node->next=NULL;
	   liste->size++;
	 }
       else
	 {
	   node->next=tmp;
	   tmp->previous=node;
	   node->previous=NULL;
	   liste->size++;
	   liste->first=node;
	   
	 }
     }
   else
     if (liste->size>1)
       {
	 while(i<liste->size-1)
	   {
	     if (tmp->inf_bound_value > value)
	       break;
	       
	     tmp=tmp->next;
	     i++;
	   }
	 if (i==liste->size-1)
	   {
	     if (tmp->inf_bound_value > value)
	       {
		 node->previous=tmp->previous;
		 node->next=tmp;
		 tmp->previous->next = node;
		 tmp->previous=node;
		 liste->size++;
	       }
	     else
	       {
		 tmp->next=node;
		 node->previous=tmp;
		 node->next=NULL;
		 liste->size++;
	       }
	   }
	 else
	   {
	     
	     if (i==0)
	       {
		node->next=tmp;
		node->previous=NULL;
		tmp->previous=node;
		liste->size++;
		liste->first=node;
		
	       }
	     else
	       {
		 node->next=tmp;
		 node->previous=tmp->previous;
		 tmp->previous->next=node;
		 tmp->previous=node;
		 liste->size++;
	       }
	   }
       }
     else
       if(liste->size==0)
	 {
	   node->next = NULL;
	   node->previous= NULL;
	   liste->first = node;
	 
	   liste->size++;
	 }


   
}

void suppress_first_liste(Liste liste)

{
  if (liste->size == 0)
    {
      printf ("\n suppress :Could not supress node : NULL liste.\n");
    }
  

  
  if (liste->size == 1)
    {
      liste->first=NULL;
      liste->size=0;
    }
  else
    {
      Node to_delete = liste->first;
      liste->first = liste->first->next;
      liste->first->previous=NULL;
      liste->size--;
      
    }
}

void clean_liste(Liste liste, double value)

{
  Node tmp=liste->first;
 
  int i=0;

  if (liste->size ==1)
    {
      if (tmp->inf_bound_value > value)
	{
	  liste->first=NULL;
	  liste->size=0;
	}
    }
  if (liste->size>1)
    {
      while (i< liste->size)
	{
	  if (tmp->inf_bound_value > value)
	    if (i==0)
	      {
		liste->first=NULL;
		liste->size=0;
	      }
	    else
	    {
	      tmp->previous->next=NULL;
	      liste->size=i;
	      break;
	    }
	  tmp=tmp->next;
	  i++;
	}
     
    }


 }


void print_liste(Liste liste)

{
  if (liste->size == 0)
    {
      printf ("\n print_liste :Could not print liste : NULL liste.\n");
    }
  else
    {
      Node current=liste->first;
      int i =0;
      
      while (i< liste->size)
	{
	  printf("\nliste elem %d \n",i);
	  printf("valeur borne : %lf\n",current->inf_bound_value);
	  current=current->next;
	  i++;
	}
    }
}
