BOXLIB_HOME ?= ../../..

TOP = $(BOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
DEBUG	      = FALSE
DIM	      = 2
COMP          = g++
FCOMP         = gfortran
USE_MPI       = TRUE
#
# Base name of the executable.
#
#EBASE = Add
#EBASE = ComputeAmrDataNorms
#EBASE = DiffFab
#EBASE = DiffSameDomainRefined
EBASE = DiffSameDomainRefinedStag
#EBASE = DiffSameGrid2
#EBASE = DiffSameGrid
#EBASE = DiffSameGridRefined
#EBASE = DiffUniform
#EBASE = PltFileNormB
#EBASE = PltFileNorm
#EBASE = PltFileScalConvRate

CEXE_sources += $(EBASE).cpp

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DBL_NOLINEVALUES -DBL_PARALLEL_IO

include Make.package
include ../Make.package
include $(BOXLIB_HOME)/Src/Base/Make.package
include $(BOXLIB_HOME)/Src/Boundary/Make.package
#include $(BOXLIB_HOME)/Src/Amr/Make.package
include $(BOXLIB_HOME)/Src/Extern/amrdata/Make.package

INCLUDE_LOCATIONS += .
INCLUDE_LOCATIONS += ..
vpathdir += ..

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base
vpathdir += $(BOXLIB_HOME)/Src/Base

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Boundary
vpathdir += $(BOXLIB_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Amr
vpathdir += $(BOXLIB_HOME)/Src/Amr

INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata
vpathdir += $(BOXLIB_HOME)/Src/Extern/amrdata

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

#include $(TOP)/Tools/GNUMake/Make.defs ./Make.package

#ifeq ($(USE_ARRAYVIEW),TRUE)
#  DEFINES += -DBL_USE_ARRAYVIEW
#  DEFINES += -DBL_ARRAYVIEW_TAGBOX
#endif





