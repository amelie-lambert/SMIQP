# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpscip.la
########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libsmiqpscip_la_SOURCES = boxqp_scip.c local_sol_scip.c sbb_miqcr_scip.c miqcr_scip.c sbb_boxqp_scip.c 

libsmiqpscip_la_DEPENDENCIES = $(libsmiqpscip_la_LIBADD)

# This is for libtool
libsmiqpscip_la_LDFLAGS = $(LT_LDFLAGS)
libsmiqpscip_la_LDFLAGS += -lpthread -llapack -lblas -lgfortran -lm -L$(SCIP_LFLAGS)  -lscip -lobjscip -llpispx2 -ltpinone -lnlpi.cppad -lsoplex.linux.x86_64.gnu.opt -lgmp -lreadline -lz -ldl -lgomp
libsmiqpscip_la_LIBADD = ../../Declarations/libsmiqpdecl.la 
libsmiqpscip_la_CPPFLAGS = -I../ -I../../Declarations -I../ -I$(SCIP_CFLAGS) 
#fi


########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = ../Declarations
includecoindir += ../SdpSolver
includecoin_HEADERS = boxqp_scip.h local_sol_scip.h sbb_miqcr_scip.h miqcr_scip.h sbb_boxqp_scip.h 



