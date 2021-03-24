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



#ifndef LISTE_H
#define LISTE_H

#include "quad_prog.h"

typedef struct node * Node;

struct node
{
  double inf_bound_value;
  double *sol_x;
  MIQCP_BAB bab;
  Node next;
  Node previous;
  
};

typedef struct liste * Liste;
struct liste
{
  Node first;
  int size;
};


Liste init_liste();
void insert_liste(Liste liste, double value, double * sol, struct miqcp_bab bab, int nb_col);
void suppress_first_liste(Liste liste);
void print_liste(Liste liste);
void clean_liste(Liste liste, double value);
#endif
