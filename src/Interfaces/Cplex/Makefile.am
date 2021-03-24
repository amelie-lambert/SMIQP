# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpcplex.la
LAPACKDIR = /usr/include
########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libsmiqpcplex_la_SOURCES = boxqp_cplex.c local_sol_cplex.c sbb_miqcr_cplex.c miqcr_cplex.c obbt.c sbb_boxqp_cplex.c

libsmiqpcplex_la_DEPENDENCIES = $(libsmiqpcplex_la_LIBADD)

# This is for libtool
libsmiqpcplex_la_LDFLAGS = $(LT_LDFLAGS)
libsmiqpcplex_la_LDFLAGS += -lpthread -llapack -lblas -lgfortran -lm -ldl -L$(CPLEX_LFLAGS) -lcplex
libsmiqpcplex_la_LIBADD = ../../Declarations/libsmiqpdecl.la ../Ipopt/libsmiqpipopt.la
libsmiqpcplex_la_LIBADD += ../Ipopt/libsmiqpipopt.la
if HAVE_SCIP
libsmiqpcplex_la_LIBADD +=../Scip/libsmiqpscip.la
endif

libsmiqpcplex_la_CPPFLAGS = -I$(CPLEX_CFLAGS) -I../ -I../../Declarations -I../Ipopt/

if HAVE_SCIP
libsmiqpcplex_la_CPPFLAGS +=-I$(SCIP_CFLAGS) -I../Scip/  -DHAVE_SCIP
endif


########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = ../../Declarations ../Ipopt
includecoin_HEADERS = boxqp_cplex.h local_sol_cplex.h sbb_miqcr_cplex.h miqcr_cplex.h obbt.h sbb_boxqp_cplex.h



