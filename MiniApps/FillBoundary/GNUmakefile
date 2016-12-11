BOXLIB_HOME = ../../../BoxLib

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_IPM      = FALSE
PROFILE      = TRUE
COMM_PROFILE = TRUE
USE_DUMPI    = FALSE
COMP         = g++
FCOMP        = gfortran
#COMP         = Intel
#FCOMP        = Intel
DIM          = 3
PRECISION    = DOUBLE
EBASE        = fbtest

CEXE_sources += MultiFabFillBoundaryTest.cpp

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs
include $(BOXLIB_HOME)/Src/Boundary/Make.package
include $(BOXLIB_HOME)/Src/Base/Make.package

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
