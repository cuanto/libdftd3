################################################################################
#
# Program makefile
#
# Needs following variables
#
# FC  Fortran compiler
# FCFLAGS  Fortran compiler flags
# LN  Linker
# LNFLAGS  Linker flags
#
################################################################################
.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: clean distclean

LIBDIR = ../lib

all: testapi testapi2

OBJS_TESTAPI = testapi.o
OBJS_TESTAPI2 = testapi2.o

testapi: $(OBJS_TESTAPI)
	$(LN) $(LNFLAGS) -o $@ $^ $(LIBDIR)/libdftd3.a

testapi2: $(OBJS_TESTAPI2)
	$(LN) $(LNFLAGS) -o $@ $^ $(LIBDIR)/libdftd3.a

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f testapi testapi2

%.o: %.f90
	$(FC) $(FCFLAGS) -I$(LIBDIR) -c $< -o $@


# Dependencies
testapi.o:
testapi2.o:
