BOXLIB_HOME = ../../../BoxLib

PRECISION = DOUBLE
DEBUG     = FALSE
DEBUG     = TRUE
DIM       = 3
COMP      = g++
FCOMP     = gfortran

USE_MPI=FALSE
USE_MPI=TRUE

USE_CXX11 = TRUE

#PROFILE       = TRUE
#TRACE_PROFILE = TRUE
#COMM_PROFILE  = FALSE

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

EBASE = iotest

HERE = .

include $(BOXLIB_HOME)/Src/Base/Make.package

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Base

CEXE_sources += IOTestDriver.cpp IOTest.cpp

vpath %.H   $(HERE) $(BOXLIB_HOME)/Src/Base
vpath %.cpp $(HERE) $(BOXLIB_HOME)/Src/Base
vpath %.F   $(HERE) $(BOXLIB_HOME)/Src/Base
vpath %.f   $(HERE) $(BOXLIB_HOME)/Src/Base
vpath %.f90 $(HERE) $(BOXLIB_HOME)/Src/Base

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

