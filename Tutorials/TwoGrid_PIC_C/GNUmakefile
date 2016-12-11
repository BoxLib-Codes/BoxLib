BOXLIB_HOME ?= ../../

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE
USE_MPI   = TRUE

USE_OMP   = FALSE

include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(BOXLIB_HOME)/Src/Base/Make.package
include $(BOXLIB_HOME)/Src/Particle/Make.package
include $(BOXLIB_HOME)/Src/Boundary/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_to_F_MG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/FParallelMG.mak
include $(BOXLIB_HOME)/Src/F_BaseLib/FParallelMG.mak

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
