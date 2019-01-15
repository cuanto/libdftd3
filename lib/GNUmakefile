################################################################################
#
# Library makefile
#
# Needs following variables
#
# FC  Fortran compiler
# FCFLAGS  Fortran compiler flags
# LN  Linker
# LNFLAGS  Linker flags
# SRCDIR  Folder, where source files can be found
#
################################################################################
.SUFFIX:
.SUFFIX: .f90 .o

LIBD3 = libdftd3.a
LIBD3SO = libdftd3.so

all: $(LIBD3) $(LIBD3SO)

vpath % $(SRCDIR)

.PHONY: clean distclean


OBJS  = api.o common.o core.o pars.o sizes.o wrapper.o

$(LIBD3): $(OBJS)
	ar r $@ $^

$(LIBD3SO): $(OBJS)
	$(FC) -shared -o $(LIBD3SO) $(OBJS) 

clean:
	rm -f *.o

distclean: clean
	rm -f *.mod $(LIBD3) $(LIBD3SO)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Dependencies
api.o: common.o core.o sizes.o
common.o:
core.o: common.o pars.o sizes.o
pars.o: common.o sizes.o
sizes.o:
