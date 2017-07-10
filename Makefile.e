# -*- Makefile -*-

##############################################################################

### F90      = 
### F90FLAGS = 
### INCLUDES = 
### LIBS     = 

LIBS0     = $(LIBS)
INCLUDES0 = $(INCLUDES)

LLIBS     = 
LINCLUDES = 

##############################################################################

include main.mk

install: all
	cp -pf $(PROG) ../../../bin/.

# list of dependencies (via USE statements)
include depend.mk
