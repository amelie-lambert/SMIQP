# Copyright (C) 2019  Amelie LAMBERT
# This file is part of the src of "Smiqp".
# Author: Amélie Lambert Cedric-CNAM-Paris 2019-10-23


AUTOMAKE_OPTIONS = foreign

CFLAGS=-g
CPPFLAGS=-g
SUBDIRS = Csdp
SUBDIRS += ConicBundle
#if COIN_HAS_IPOPT
SUBDIRS += Ipopt
#endif
SUBDIRS += src/Declarations
SUBDIRS += src/SdpSolver
SUBDIRS += src/Interfaces
#if COIN_HAS_IPOPT
SUBDIRS += src/Interfaces/Ipopt
#endif

if HAVE_SCIP
SUBDIRS += src/Interfaces/Scip
endif
if HAVE_CPLEX
SUBDIRS += src/Interfaces/Cplex
endif
SUBDIRS += src/Alg


EXTRA_DIST = autogen.sh

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool



ACLOCAL_AMFLAGS = '-I m4'


########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into

#doc_DATA = README.md AUTHORS LICENSE
