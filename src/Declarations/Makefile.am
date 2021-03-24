# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpdecl.la

########################################################################
#                          smiqp library                              #
########################################################################

# List all source files, including headers
libsmiqpdecl_la_SOURCES = liste.c quad_prog.c quicksort.c utilities.c 

# This is for libtool
libsmiqpdecl_la_LDFLAGS = $(LT_LDFLAGS)
libsmiqpdecl_la_LDFLAGS += -llapack -lblas -lgfortran -lm 

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = liste.h quad_prog.h quicksort.h




