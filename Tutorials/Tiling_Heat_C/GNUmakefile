# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../..

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = TRUE
PROFILE   = FALSE
COMP      = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = Heat

include ./Make.package
include $(BOXLIB_HOME)/Tools/GNUMake/Make.defs
include $(BOXLIB_HOME)/Src/Base/Make.package

include $(BOXLIB_HOME)/Tools/GNUMake/Make.rules
