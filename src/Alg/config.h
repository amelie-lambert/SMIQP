/* src/Alg/config.h.  Generated from config.h.in by configure.  */
/* src/Alg/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if ASL is available. */
/* #undef COIN_HAS_ASL */

/* Define to 1 if Blas is available. */
#define COIN_HAS_BLAS 1

/* Define to 1 if the LAPACK package is available */
#define COIN_HAS_LAPACK 1

/* Define to a macro mangling the given C identifier (in lower and upper
   case). */
#define COIN_LAPACK_FUNC(name,NAME) name ## _

/* As COIN_LAPACK_FUNC, but for C identifiers containing underscores. */
#define COIN_LAPACK_FUNC_(name,NAME) name ## _

/* Define to the debug sanity check level (0 is no test) */
/* #undef COIN_MIQCR_CHECKLEVEL */

/* Define to the debug verbosity level (0 is no output) */
/* #undef COIN_MIQCR_VERBOSITY */

/* Define to 1 if your C++ compiler doesn't accept -c and -o together. */
/* #undef CXX_NO_MINUS_C_MINUS_O */

/* Define to the C type corresponding to Fortran INTEGER */
#define FORTRAN_INTEGER_TYPE int

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Version number of project */
/* #undef MIQCR_VERSION */

/* Major version number of project. */
/* #undef MIQCR_VERSION_MAJOR */

/* Minor version number of project. */
/* #undef MIQCR_VERSION_MINOR */

/* Release version number of project. */
/* #undef MIQCR_VERSION_RELEASE */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://cedric.cnam.fr/~lamberta/smiqp"

/* Define to the full name of this package. */
#define PACKAGE_NAME "smiqp"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "smiqp 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "smiqp"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1
