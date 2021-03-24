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

#ifndef SOLVER_SDP_H
#define SOLVER_SDP_H

void compute_alpha_beta_miqcr() ;
int gen_subg(void *function_key,double *X,int n_indices,int *variable_indices,double *new_subgradient_values);
void run_conic_bundle(double * obj_val);

#endif
