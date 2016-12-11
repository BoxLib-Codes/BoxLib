BOXLIB_HOME ?= ../../../BoxLib

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
USE_IPM      = FALSE
PROFILE      = TRUE
TRACE_PROFILE = TRUE
COMM_PROFILE = TRUE
USE_DUMPI    = FALSE
COMP         = g++
DIM          = 3
PRECISION    = DOUBLE
EBASE        = mkdirtest

CEXE_sources += MKDir.cpp

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs
include $(BOXLIB_HOME)/Src/Base/Make.package

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
