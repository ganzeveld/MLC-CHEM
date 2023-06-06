##############################################################################

### choose Fortran90 compiler and options

SYSTEM := $(shell uname)
ifeq ($(SYSTEM),OSF1)
  ### f90 on alphas
  F90       = f90
  F90FLAGS  = -O -cpp
  F90R8     = -r8
else
  ifeq ($(SYSTEM),Linux)
    ### lf95 = Lahey
    F90 = lf95
    ### --dbl = double precision (necessary because some variables use _dp)
    ### -Cpp  = run C-preprocessor before compiling
    ###         (necessary because of some compiler directives in the code)
    ### --pca = protect constant argument (necessary for unknown reason)
    #F90FLAGS = -Cpp --pca
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0
    #F90FLAGS = -Cpp --chk a,e,s   --pca --ap -O0 -g --trap --verbose
    F90FLAGS = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --verbose
    F90R8    = --dbl
  else
	# IBM-Darwin / MAC OS X (IBM XL F/C Compiler)
	ifeq ($(SYSTEM),Darwin)
	  F90      = xlf95_r
	  F90FLAGS = -O2 -pg -qarch=auto -qtune=auto -qextname -qsuppress=1518-061:1518-128 -qstrict -qfixed=100 -qMAXMEM=-1 -qsuffix=f=f90 -qsuffix=cpp=f90 -qzerosize -WF,-D__ibm__ -d -WF,-qlanglvl=classic -qlanglvl=95pure -qfree=f90
	  F90R8    = -qrealsize=8
	else
		ERROR
	endif
  endif
endif

# activate following lines for intel compiler (8.0.039) to overwrite
# previous settings
F90       = ifort
#F90FLAGS  = -cpp -O0
# more extensive compiler options to do full check, activate by removing # for testing array boundaries etc.
F90FLAGS  = -cpp -i4 -r8 -O3 -align dcommons -check bounds -extend_source
F90R8     = -autodouble

# activate following lines for g95 to overwrite previous settings
#F90       = g95
#F90FLAGS  = -cpp -O0
#F90R8     = -r8

##############################################################################

### The above block defines the variables F90, F90FLAGS, INCLUDES, and
### LIBS. The if...else...endif constructs try to find suitable values
### for different architectures and machines. To implement the mecca
### boxmodel on a new machine, you can simply overwrite the above
### definitions by activating the following block and entering
### appropriate values:
### F90       = 
### F90FLAGS  = 
### F90R8     =

##############################################################################

# targets
include main.mk

# list of dependencies (via USE statements)
include depend.mk
