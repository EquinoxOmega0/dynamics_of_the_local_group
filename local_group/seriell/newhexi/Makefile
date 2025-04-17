.KEEP_STATE:
.SUFFIXES: .f90 .F90

#=#############################################
#            sources (.f90, .F90)
#=#############################################

# source-directory
SRCDIR = src

# source-files
SOURCE_FIRST = \
nrtype.f90

# source-files
SOURCE = \
nr.f90 nrutil.f90 recipes.f90 newhexi.f90

#=#############################################
#             general library definitions
#=#############################################

#--------------------- 
# FORTRAN-definitions
#--------------------- 

#=#############################################
#             LINUX definitions
#=#############################################

#-----------------------
#         ifort
#-----------------------

FORTRAN = ifort
FFLAGS  = -O2

#--------------------- 
#  libraries
#--------------------- 
 
LIBS =
 
#--------------------- 
#  object files
#--------------------- 
 
OBJECTS = $(SOURCE_FIRST:.f90=.o) $(SOURCE:.f90=.o)
 
#================================================================
#                          targets
#================================================================

default:	$(OBJECTS)
	$(FORTRAN) $(OBJECTS) $(FFLAGS) -o newhexi.out $(LIBS)

clean:
	rm -f *.o *.mod .m* .n* *.lst ../.make.state ../.nse_depinfo

.f90.o:
	$(FORTRAN) $(FFLAGS) -c $<

.F90.o:
	$(FORTRAN) $(FFLAGS) -c $<
