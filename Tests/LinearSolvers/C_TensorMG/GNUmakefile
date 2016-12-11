BOXLIB_HOME := ../../..

PRECISION = DOUBLE

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 2
DIM	= 3

COMP    = gnu

USE_MPI = FALSE

EBASE = main

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs

CXXFLAGS +=
CXXOPTF  +=
CXXDEBF  +=

include ./Make.package

include $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/Base/Make.package
include $(BOXLIB_HOME)/Src/Boundary/Make.package

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_TensorMG

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Boundary
vpathdir          += $(BOXLIB_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base
vpathdir          += $(BOXLIB_HOME)/Src/Base

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Amr
vpathdir          += $(BOXLIB_HOME)/Src/Amr

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG
vpathdir          += $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG

CEXE_sources += BCRec.cpp 
CEXE_headers += BCRec.H BC_TYPES.H

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules


