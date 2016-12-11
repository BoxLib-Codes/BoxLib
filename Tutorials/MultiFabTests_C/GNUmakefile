# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../../../BoxLib

DEBUG     = FALSE
DEBUG     = TRUE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = TRUE
COMM_PROFILE  = TRUE
TRACE_PROFILE = TRUE
COMP      = g++
#COMP      = Intel
DIM       = 3
PRECISION = DOUBLE
EBASE     = mfreadwrite
EBASE     = gridmovetest

include ./Make.package
include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs
include $(BOXLIB_HOME)/Src/Boundary/Make.package
include $(BOXLIB_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Boundary
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Amr

vpathdir += $(BOXLIB_HOME)/Src/Boundary
vpathdir += $(BOXLIB_HOME)/Src/Base
vpathdir += $(BOXLIB_HOME)/Src/Amr

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
