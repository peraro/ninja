#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Ninja], [1.0.0], [tiziano.peraro@ed.ac.uk])
AC_CONFIG_SRCDIR([src/ninja_implem.cxx])
AC_CONFIG_HEADERS([config.h src/ninja/ninja_config.h])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/m4_ax_cxx_compile_stdcxx_11.m4])
AM_INIT_AUTOMAKE
LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_FC

# From now on is c++ stuff
AC_LANG_PUSH([C++])

# AC_OPENMP
ACX_PTHREAD

# Check C++11 compiler support
AX_CXX_COMPILE_STDCXX_11
AX_CXX_COMPILE_STDCXX_11_INIT_LIST
AX_CXX_COMPILE_STDCXX_11_RANDOM

# this is only defined in config.h
AC_DEFINE([NINJA_CONFIG_H_INTERNAL],
          [1],
          [Define to 1 if this an internal is a source file.])

# Choose whether the interface to the OneLoop library should be
# compiled
AC_ARG_WITH([avholo],
    AS_HELP_STRING([--with-avholo@<:@=FLAGS@:>@],
    [include an interface with the OneLoop integral library,\
    specifying the corresponding flags for dynamic linking \
    (if the Fortran module 'avh_olo' is not in a standard path, \
    you should add its directory to the FCINCLUDE variable)]),
    [],
    [with_avholo=no])
AS_IF([test "x$with_avholo" == xyes],
       [AC_SUBST([AVHFLAGS], ["-lavh_olo"])
        AC_DEFINE([NINJA_USE_ONELOOP], [1],
        [Define to 1 if you have OneLoop.])
        [avholo=true]
        ],
       [test "x$with_avholo" != xno],
       [AC_SUBST([AVHFLAGS], ["$with_avholo"])
        AC_DEFINE([NINJA_USE_ONELOOP], [1],
        [Define to 1 if you have OneLoop.])
        [avholo=true]
        ],
        [avholo=false])
AM_CONDITIONAL(AVHOLO, test "x$avholo" == xtrue)

# Choose whether an internal integral cache for the OneLoop library
# should be used
AC_ARG_ENABLE([avholo_cache],
    AS_HELP_STRING([--disable-avholo_cache],
    [do not include a cache of Master Integrals for OneLoop]))
AS_IF([test "x$enable_avholo_cache" != xno],
      [AC_DEFINE([NINJA_USE_ONELOOP_WITH_CACHE], [1],
      [Define to 1 if a cache of Master Integrals is implemented for OneLoop.])
      ],
      [])

# Choose whether the interface to the LoopTools library should be
# compiled
AC_ARG_WITH([looptools],
   AS_HELP_STRING([--with-looptools@<:@=FLAGS@:>@],
                  [include an interface with the LoopTools library, \
                  specifying the corresponding flags for static linking. \
                  (if the header file cloopttols.h is not in a standard \
                  path, you should add its directory to the \
                  CPPFLAGS variable)]),
                  [],
                  [with_looptools=no])
AS_IF([test "x$with_looptools" == xyes],
       [AC_SUBST([LTFLAGS], ["-looptools"])
        AC_DEFINE([NINJA_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true]
        ],
       [test "x$with_looptools" != xno],
       [AC_SUBST([LTFLAGS], ["$with_looptools"])
        AC_DEFINE([NINJA_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true]
        ],
        [looptools=false])
AM_CONDITIONAL(LOOPTOOLS, test "x$looptools" == xtrue)
AS_IF([test "x$looptools" == xtrue],
    [[case $host_os in
      darwin*)
         exninjastatic=true
    esac]])
AM_CONDITIONAL(EXNINJASTATIC, test "x$exninjastatic" == xtrue)

# Always enable massless specialization
AC_DEFINE([NINJA_MASSLESS], [1],
          [Define to 1 massless specialization is enabled.])

# Choose whether to include higher-rank support
AC_ARG_ENABLE([higher_rank],
    AS_HELP_STRING([--enable-higher_rank],
                   [enable support for higher-rank numerators]))
AS_IF([test "x$enable_higher_rank" != xyes],[],
      [AC_DEFINE([NINJA_X1RANK], [1],
      [Define to 1 if higher-rank support is enabled.])])

# Choose whether to compile in quadruple precision
AC_ARG_WITH([quadruple],
    AS_HELP_STRING([--with-quadruple@<:@=FLAGS@:>@],
    [compile the library in quadruple precision (requires GCC libquadmath), \
    specifying the corresponding flags for the linker]),
    [],
    [with_quadruple=no])
AS_IF([test "x$with_quadruple" == xyes],
       AC_SUBST([QUADRUPLEFLAGS], ["-lquadmath"])
       AC_DEFINE([NINJA_QUADRUPLE], [1],
       [Define to 1 if Ninja is compiled in quadruple precision.])
       [quadruple=true],
       [test "x$with_quadruple" != xno],
       AC_SUBST([QUADRUPLEFLAGS], ["$with_quadruple"])
       AC_DEFINE([NINJA_QUADRUPLE], [1],
       [Define to 1 if Ninja is compiled in quadruple precision.])
       [quadruple=true],
       [quadruple=false])
AM_CONDITIONAL(QUADRUPLE, test "x$quadruple" == xtrue)

# Choose whether compiling interfaces: gosam
AC_ARG_ENABLE([gosam],
    AS_HELP_STRING([--disable-gosam], [do not include GoSam interface]))
AS_IF([test "x$enable_gosam" != xno], [gosam=true],
 [gosam=false]
  )
AM_CONDITIONAL(GOSAM, test "x$gosam" == xtrue)
AM_CONDITIONAL(NINJAGOSAMCONTRIB, false)

AC_ARG_VAR([FCINCLUDE],[additional Fortran include flags, e.g. \
                        -I<include dir> if you have fortran modules \
                        installed in a nonstandard directory <include dir>])

# Checks for typedefs, structures, and compiler characteristics.
# AC_CHECK_HEADER_STDBOOL
# AC_C_INLINE
# AC_TYPE_SIZE_T

# Checks for library functions.
# AC_CHECK_FUNCS([floor pow rint sqrt])


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/ninja/Makefile
				 examples/Makefile])
AC_OUTPUT
