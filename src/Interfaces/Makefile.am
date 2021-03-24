# This file is distributed under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory

noinst_LTLIBRARIES = libsmiqpinterface.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libsmiqpinterface_la_SOURCES = in_out.c parameters.c

libsmiqpinterface_la_DEPENDENCIES = $(libsmiqpinterface_la_LIBADD)
libsmiqpinterface_la_LDFLAGS = $(LT_LDFLAGS)
# This is for libtool
libsmiqpinterface_la_LDFLAGS += -lm
libsmiqpinterface_la_LIBADD = ../Declarations/libsmiqpdecl.la
libsmiqpinterface_la_CPPFLAGS = -I../Declarations

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = ../Declarations
includecoin_HEADERS = in_out.h parameters.h



