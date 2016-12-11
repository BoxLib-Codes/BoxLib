
DIM          = 3

COMP         = gnu

DEBUG        = FALSE

USE_MPI      = TRUE
USE_OMP      = FALSE
USE_UPCXX    = TRUE

BOXLIB_HOME = ../..

EBASE = main

include ./Make.package

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs

include $(BOXLIB_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base

vpathdir += $(BOXLIB_HOME)/Src/Base

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
