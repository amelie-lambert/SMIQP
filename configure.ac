# This file is distributed under under the  GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

## $Id: configure.ac 2595 2020-01-17 22:10:23Z  $

# Authors:  Amelie Lambert    2020-01-17

#############################################################################
#                       Names and other basic things                        #
#############################################################################


AC_CONFIG_MACRO_DIRS([BuildTools/])

AC_PREREQ(2.59)
AC_INIT([smiqp],[1.0],[http://cedric.cnam.fr/~lamberta/smiqp])

AC_COPYRIGHT(
[Copyright 2019. Amélie Lambert. All Rights Reserved.
This file is part of the open source package SMIQP which is distributed
under the Eclipse Public License.])


# List one file in the package so that the configure script can test
# whether the package is actually there

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/Alg/smiqp.c)

# Do some project-level initialization work (version numbers, ...)
AC_COIN_INITIALIZE(Smiqp)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C, C++, and Fortran compilers and appropriate compiler options.
AC_COIN_PROG_CC
AC_COIN_PROG_CXX
AC_PROG_CXXCPP 
AC_PROG_F77
AM_PROG_CC_C_O

# If there is a Fortran compiler, then setup everything to use it, including F77_FUNC
#if test -n "$F77" ; then
#  AC_F77_SETUP
#fi

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
AC_COIN_RPATH_FLAGS([$expanded_libdir])


# with a C or Fortran compiler
AC_COIN_CXXLIBS

# Doxygen
AC_COIN_DOXYGEN

# COIN_SMIQP_VERBOSITY and COIN_SMIQP_DEBUGLEVEL
AC_COIN_DEBUGLEVEL(Smiqp)

AC_ARG_WITH([cplex-cflags],
[  --with-cplex-cflags    use  cplex as subsolver],
[ use_cplex=yes ],[ use_cplex=no ])
AM_CONDITIONAL([HAVE_CPLEX], [test $use_cplex == yes])

AC_ARG_WITH([scip-cflags],
[  --with-scip-cflags    use scip as subsolver],
[ use_scip=yes],[ use_scip=no ])
AM_CONDITIONAL([HAVE_SCIP], [test $use_scip == yes])
echo "use_scip"
echo $use_scip
echo "use_cplex"
echo $use_cplex


#############################################################################
#                              COIN-OR components                           #
#############################################################################

#AC_COIN_CHK_LIBM(SmiqpLib)
AC_CONFIG_HEADER(config.h)
AC_COIN_CHK_PKG(Blas)
if test $coin_has_blas != yes; then
  AC_MSG_WARN([Warning package BLAS not found.])
fi

AC_COIN_CHK_LAPACK(SmiqpLib)
if test $coin_has_lapack != yes; then
  AC_MSG_WARN([Compiling code without LAPACK.  Certain options (e.g., quasi-Newton) will not work.])
fi

#AC_COIN_CHK_PKG(Ipopt)
#if test $coin_has_ipopt != yes; then
 #  AC_MSG_ERROR([Required package Ipopt not found.])
 #fi
 #if test $coin_has_ipopt == yes; then
AC_CONFIG_SUBDIRS([Ipopt])
#fi

AC_COIN_CHK_LIB([ConicBundle], [CB],[cb_cinterface.h])
AC_COIN_CHK_LIB([Csdp], [CSDP],[declaration.h])

AC_COIN_CHK_LIB(Cplex,[cplex],[LCPLEX=" -lcplex"],[],[INCCPLEX="ilcplex/cplex.h"])
if test $use_cplex == yes; then
   
   AC_CHECK_HEADER([ilcplex/cplex.h])
fi



 AC_COIN_CHK_LIB(Scip, [Scip],[LSCIP=" -lscip"],[],[INCSCIP1="scipdefplugins.h" INCSCIP2=scip.h"])	
if test $use_scip == yes; then
  
   AC_CHECK_HEADER([scipdefplugins.h])
   AC_CHECK_HEADER([scip.h])
 fi

echo "coin_has_scip"
echo $coin_has_scip
echo "coin_has_cplex"
echo $coin_has_cplex

#AC_COIN_CHK_LIB(Ipopt, [Ipopt],[LIPOPT=" -lipopt"],[],[INCIPOPT="IpStdCInterface.h"])
#if test $coin_has_ipopt == yes; then
#  AC_CHECK_HEADER([IpStdCInterface.h])
#fi

AC_COIN_CHK_PKG(ASL,[SmiqpAmplInterfaceLib],[coinasl],[build])


############################################################################
############################################################################
#                   Stuff that we need for C programs                    #
############################################################################
############################################################################


AC_LANG_PUSH(C)
AC_CHECK_HEADER([stdio.h])
AC_CHECK_HEADER([time.h])
AC_CHECK_HEADER([assert.h])
AC_CHECK_HEADER([stdio.h])
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([sys/stat.h])
AC_CHECK_HEADER([fcntl.h])
AC_CHECK_HEADER([unistd.h])
AC_CHECK_HEADER([sys/type.h])
AC_CHECK_HEADER([sys/wait.h])
AC_CHECK_HEADER([sys/timeb.h])
AC_CHECK_HEADER([string.h])
AC_CHECK_HEADER([lapacke.h])
AC_CHECK_HEADER([declaration.h])
AC_CHECK_HEADER([cb_cinterface.h])


######################################
# Equivalent int Fortran and C types #
######################################

# FIXME: The following test should be active, but this requires chance in
#        code to copy Index* to ipfint* arrays...

if test "$cross_compiling" = no && test "$is_bg" != yes; then
  
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int],[Define to the C type corresponding to Fortran INTEGER])
  
else
  AC_MSG_WARN([We are cross compiling, assuming Fortran `INTEGER' type corresponds to C `int' type])
  AC_DEFINE([FORTRAN_INTEGER_TYPE],[int])
fi

AC_LANG_POP(C)


#########################
# Makefile conditionals #
#########################

# The following variable collects the names of libraries that should
# be included into libsmiqp.a (relative to subdir Interfaces, where it
# is made)

#AC_SUBST(IPALLLIBS)
#IPALLLIBS="../contrib/CGPenalty/libcgpenalty.la ../Algorithm/libipoptalg.la ../Algorithm/LinearSolvers/liblinsolvers.la ../Common/libcommon.la ../LinAlg/liblinalg.la ../LinAlg/TMatrices/libtmatrices.la"


AM_INIT_AUTOMAKE

LT_PREREQ([2.2])

LT_INIT


########################################################################
##          Create Makefiles and other stuff                          ##
########################################################################

AC_CONFIG_FILES([Makefile
		 src/Declarations/Makefile
		 src/SdpSolver/Makefile
		 src/Interfaces/Makefile
		 src/Interfaces/Ipopt/Makefile	
		 src/Alg/Makefile])

if test $use_cplex == yes; then
	AC_CONFIG_FILES([src/Interfaces/Cplex/Makefile])
 fi

if test $use_scip == yes; then
AC_CONFIG_FILES([src/Interfaces/Scip/Makefile])
 fi


AC_OUTPUT

#if test $coin_has_asl = yes ; then
#  AC_CONFIG_FILES([smiqpamplinterface.pc:src/Apps/AmplSolver/smiqpamplinterface.pc.in
#                   smiqpamplinterface-uninstalled.pc:src/Apps/AmplSolver/smiqpamplinterface-uninstalled.pc.in
#                  ])
#fi

AC_CONFIG_HEADER([src/Alg/config.h])

AC_COIN_FINALIZE
