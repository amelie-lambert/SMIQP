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

#ifndef OBBT_H
#define OBBT_H


#include"quad_prog.h"


void presolve_with_obbt(int * nb_pass_var, int * nb_pass_cont,double upperbound, int nb_pass);
void presolve_with_obbt_bb(struct miqcp_bab bab,int * nb_pass_var, int * nb_pass_cont, double upperbound,int nb_pass);
void initialize_bounds_with_obbt();
#endif
