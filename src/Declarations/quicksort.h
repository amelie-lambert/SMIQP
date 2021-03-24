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

#ifndef QUICKSORT_H
#define QUICKSORT_H


struct tab_ind {
  int ind;
  double val;
};

typedef struct tab_ind * TAB_IND;

TAB_IND * create_tab_ind_from_check_violated(double *check_violated, int n);

int initialize_constraints_tab(TAB_IND ** tab,int * constraints, int * constraints_old,int nb_max_cont);
void quickSort(TAB_IND ** tab, int beg, int end);
#endif
  
