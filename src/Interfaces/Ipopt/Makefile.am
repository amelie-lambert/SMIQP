
# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpipopt.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libsmiqpipopt_la_SOURCES = local_sol_ipopt.c sbb_miqcr_ipopt.c sbb_boxqp_ipopt.c 

libsmiqpipopt_la_DEPENDENCIES = ../../../Ipopt/src/Interfaces 
libsmiqpipopt_la_LDFLAGS =  ../../../Ipopt/src/Interfaces/.libs
# This is for libtool
libsmiqpipopt_la_LDFLAGS += -lpthread -llapack -lblas -lgfortran -lm -L../../../Ipopt/src/Interfaces/.libs -lipopt
libsmiqpipopt_la_LIBADD = ../../Declarations/libsmiqpdecl.la 
libsmiqpipopt_la_CPPFLAGS = -I../../../Ipopt/src/Interfaces -I../../../Ipopt/src/Common -I../../../Ipopt/src/contrib/LinearSolverLoader -I.. -I../../Declarations 



########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = ../Declarations
includecoin_HEADERS = local_sol_ipopt.h sbb_miqcr_ipopt.h sbb_boxqp_ipopt.h 



