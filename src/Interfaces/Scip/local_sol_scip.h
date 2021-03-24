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

#ifndef LOCAL_SOL_SCIP_H
#define LOCAL_SOL_SCIP_H


/*double best_sol_adm=0;*/
int compute_local_sol_scip_init( ) ;
int compute_local_sol_scip_bb (double ** l, double ** u);


#endif
