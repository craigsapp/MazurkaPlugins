##
## Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
## Creation Date: Thu Mar  3 21:27:40 PST 2016
## Last Modified: Thu Mar  3 21:27:42 PST 2016
## Filename:      MazurkaPlugins/external/Makefile
## Syntax:        GNU Makefile
##
## Description:   This Makefile prepares the external libraries.
##
## To run this makefile, type:
##     make everything
##

.PHONY: fftw vamp

all:
	@echo
	@echo "fftw       == compile FFTW"
	@echo "vamp       == compile Vamp SDK"
	@echo "clean      == clean external libraries"
	@echo "superclean == really clean external libraries"
	@echo


everything: fftw vamp


fftw:
	(cd fftw; $(MAKE) everything)


vamp:
	(cd vamp; $(MAKE) everything)


clean:
	(cd fftw; $(MAKE) clean)
	(cd vamp; $(MAKE) clean)


superclean: clean
	(cd fftw; $(MAKE) superclean)
	(cd vamp; $(MAKE) superclean)


