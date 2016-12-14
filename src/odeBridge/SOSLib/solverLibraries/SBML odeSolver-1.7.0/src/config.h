/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Major Version of GRAPHVIZ Library */
/* #undef GRAPHVIZ_MAJOR_VERSION */

/* Minor Version of GRAPHVIZ Library */
/* #undef GRAPHVIZ_MINOR_VERSION */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

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

/* Name of package */
#define PACKAGE "SBML_odeSolver"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "{raim,xtof}@tbi.univie.ac.at"

/* Define to the full name of this package. */
#define PACKAGE_NAME "odeSolver"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "odeSolver 1.7.0beta"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "SBML_odeSolver"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.7.0beta"

/* SBML include directories */
#define SBML_CPPFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib/include/sbml"

/* SBML include directories */
#define SBML_CPPFLAGS2 "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib/include/sbml"

/* SBML lib directories */
#define SBML_LDFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib/lib"

/* SBML libs */
#define SBML_LIBS "sbml"

/* shared library extrension */
#define SHAREDLIBEXT ".so"

/* SOSLIB include directories */
#define SOSLIB_CPPFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/include"

/* SOSLIB lib directories */
#define SOSLIB_LDFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib"

/* SOSLIB libs */
#define SOSLIB_LIBS "ODES"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* SUNDIALS include directories */
#define SUNDIALS_CPPFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib/include"

/* SUNDIALS lib directories */
#define SUNDIALS_LDFLAGS "/home/chuanfuyap/NetbeansProjects/SOSlibJNIC/solverLibraries/lib/lib"

/* SUNDIALS libs */
#define SUNDIALS_LIB1 "sundials_ida"

/* SUNDIALS libs */
#define SUNDIALS_LIB2 "sundials_kinsol"

/* SUNDIALS libs */
#define SUNDIALS_LIB3 "sundials_cvodes"

/* SUNDIALS libs */
#define SUNDIALS_LIB4 "sundials_nvecserial"

/* SUNDIALS libs */
#define SUNDIALS_LIB5 "sundials_shared"

/* SUNDIALS libs */
#define SUNDIALS_LIB6 "m"

/* Define to 1 to use the GRACE Library */
#define USE_GRACE 0

/* Define to 1 to use the GRAPHVIZ Library */
#define USE_GRAPHVIZ 0

/* Define to 1 to use the SUNDIALS Library */
#define USE_SUNDIALS 1

/* Version number of package */
#define VERSION "1.7.0beta"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
