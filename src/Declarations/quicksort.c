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
#include<math.h>
#include"utilities.h"
#include <time.h>
#include <stdlib.h>
#include"quicksort.h"

extern double EPS_VIOL_CB;

TAB_IND * create_tab_ind_from_check_violated(double * check_violated, int n)
{
  int i;
  TAB_IND * tab = (TAB_IND*) malloc(n*sizeof(struct tab_ind *));
  for(i=0;i<n;i++)
    tab[i]= (TAB_IND) malloc(sizeof(struct tab_ind));

  for(i=0;i<n;i++)
    {
      tab[i]->val = check_violated[i];
      tab[i]->ind = i;
    }
  

  return tab;
}


int generate_pivot(int beg, int end){
  srand(time(NULL));
  int pivot = rand() % (end - beg)+ beg;
  return pivot;
}
 


void change(TAB_IND ** tab, int i, int j)
{
  TAB_IND temp;
  temp = tab[0][i];
  tab[0][i] = tab[0][j];
  tab[0][j] = temp;
}

int partition(TAB_IND ** tab, int beg, int end, int pivot) {
  int i, j;
  
  change(tab,pivot,end);
  j=beg; 

  for(i=beg;i<end;i++)
    {
      if (tab[0][i]->val - tab[0][end]->val < -EPS_VIOL_CB){
	change(tab,i,j);
	j++;

      }
    }
  change(tab,end,j);
  
  return j;
}

void quickSort(TAB_IND ** tab, int beg, int end)
{
  int pivot,i;
   if( beg < end ) 
   {
     pivot = generate_pivot(beg, end);
     pivot = partition( tab, beg, end,pivot);
     quickSort( tab, beg, pivot-1);
     quickSort( tab, pivot+1, end);
   }


   
  
}

int initialize_constraints_tab(TAB_IND ** tab,int * constraints, int * constraints_old, int nb_max_cont)
{
  int i;
  
  
  for (i=0; i< nb_max_cont; i++)
    {
      if (tab[0][i]->val > EPS_VIOL_CB)
	break;

      
      constraints[i]=tab[0][i]->ind+1;
      constraints_old[i]=tab[0][i]->ind+1;
    }

  
  return i;
}
