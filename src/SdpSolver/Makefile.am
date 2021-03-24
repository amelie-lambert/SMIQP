# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpsdp.la

########################################################################
#                          smiqp library                      #         #
# ########################################################################

# List all source files, including headers
libsmiqpsdp_la_SOURCES = solver_sdp.c solver_sdp_mixed.c

libsmiqpsdp_la_DEPENDENCIES = $(libsmiqpsdp_la_LIBADD)
libsmiqpsdp_la_LIBADD = ../Declarations/libsmiqpdecl.la

libsmiqpsdp_la_LDFLAGS = $(LT_LDFLAGS)
# This is for libtool
libsmiqpsdp_la_LDFLAGS += -L../../ConicBundle/lib -lcb -L../../Csdp/lib -lsdp -llapack -lblas -lgfortran -lm 
libsmiqpsdp_la_CPPFLAGS = -I../../ConicBundle/include -I../../Csdp/include -I../Declarations/ -g -m64 -fPIC -D_GLIBC_DEBUG -DNOSHORTS -W

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = ../Interfaces
includecoindir += ../Declarations
includecoin_HEADERS = solver_sdp.h solver_sdp_mixed.h



