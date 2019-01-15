.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: dftd3 lib testapi testapi2

all: lib dftd3 testapi testapi2

include make.arch

lib:
	$(MAKE) -C lib CC="$(CC)" FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
            LNFLAGS="$(LNFLAGS)" SRCDIR="."

testapi: lib
	$(MAKE) -C test FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	    LNFLAGS="$(LNFLAGS)" testapi

testapi2: lib
	$(MAKE) -C test FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	    LNFLAGS="$(LNFLAGS)" testapi2


.PHONY: clean distclean
clean:
	$(MAKE) -C lib clean
	$(MAKE) -C prg clean
	$(MAKE) -C test clean

distclean:
	$(MAKE) -C lib distclean
	$(MAKE) -C prg distclean
	$(MAKE) -C test distclean
