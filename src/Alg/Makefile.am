# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17


########################################################################x
#                                smiqp                                 #
########################################################################






#include_HEADERS = ../Interfaces/parameters.h ../Declarations/utilities.h ../Declarations/liste.h ../Declarations/quicksort.h ../Interfaces/Cplex/obbt.h ../Interfaces/in_out.h ../Interfaces/Cplex/local_sol_cplex.h ../Interfaces/Scip/local_sol_scip.h ../Interfaces/Ipopt/local_sol_ipopt.h ../Interfaces/Cplex/miqcr_cplex.h ../Interfaces/Cplex/sbb_miqcr_cplex.h ../Interfaces/Ipopt/sbb_miqcr_ipopt.h ../Interfaces/Scip/miqcr_scip.h ../Interfaces/Scip/sbb_miqcr_scip.h ../Interfaces/Scip/sbb_boxqp_scip.h ../Interfaces/Ipopt/sbb_boxqp_ipopt.h ../Interfaces/Cplex/sbb_boxqp_cplex.h ../Interfaces/Cplex/boxqp_cplex.h ../Interfaces/Scip/boxqp_scip.h ../SdpSolver/solver_sdp.h ../SdpSolver/solver_sdp_mixed.h ../Declarations/quad_prog.h


#SUBDIRS = ../Declarations ../Interfaces/ ../Interfaces/Scip/ ../Interfaces/Ipopt ../Interfaces/Cplex/ ../include

bin_PROGRAMS = smiqp
CCLD = g++
smiqp_SOURCES = smiqp.c
#../Interfaces/parameters.h ../Declarations/utilities.h ../Declarations/liste.h ../Declarations/quicksort.h ../Interfaces/Cplex/obbt.h ../Interfaces/in_out.h ../Interfaces/Cplex/local_sol_cplex.h ../Interfaces/Scip/local_sol_scip.h ../Interfaces/Ipopt/local_sol_ipopt.h ../Interfaces/Cplex/miqcr_cplex.h ../Interfaces/Cplex/sbb_miqcr_cplex.h ../Interfaces/Ipopt/sbb_miqcr_ipopt.h ../Interfaces/Scip/miqcr_scip.h ../Interfaces/Scip/sbb_miqcr_scip.h ../Interfaces/Scip/sbb_boxqp_scip.h ../Interfaces/Ipopt/sbb_boxqp_ipopt.h ../Interfaces/Cplex/sbb_boxqp_cplex.h ../Interfaces/Cplex/boxqp_cplex.h ../Interfaces/Scip/boxqp_scip.h ../SdpSolver/solver_sdp.h ../SdpSolver/solver_sdp_mixed.h ../Declarations/quad_prog.h

smiqp_LDADD =  ../Declarations/libsmiqpdecl.la ../SdpSolver/libsmiqpsdp.la ../Interfaces/libsmiqpinterface.la ../Interfaces/Ipopt/libsmiqpipopt.la -lm -ldl -lpthread -llapack -lblas -lgfortran 

smiqp_CPPFLAGS = -I../Declarations -I../SdpSolver -I../Interfaces -I../Interfaces/Ipopt  -I../include -lm -ldl -lpthread -llapack -lblas -lgfortran  -llpispx2 -ltpinone -lgmp -lreadline -lz  -lgomp  -llpispx2 -ltpinone -lgmp -lreadline -lz  -lgomp



if HAVE_CPLEX
smiqp_CPPFLAGS +=-I../Interfaces/Cplex -DHAVE_CPLEX -L$(CPLEX_LFLAGS) -lcplex
smiqp_LDADD +=../Interfaces/Cplex/libsmiqpcplex.la
endif

if HAVE_SCIP
smiqp_CPPFLAGS += -I../Interfaces/Scip -DHAVE_SCIP -L$(SCIP_LFLAGS)  -lscip -lobjscip -lnlpi.cppad -lsoplex.linux.x86_64.gnu.opt
smiqp_LDADD += ../Interfaces/Scip/libsmiqpscip.la

endif

smiqp_LDADD +=  -lm -ldl -lpthread -llapack -lblas -lgfortran 

