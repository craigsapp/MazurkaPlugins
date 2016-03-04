##
## Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
## Creation Date: Thu Mar  3 20:57:41 PST 2016
## Last Modified: Thu Mar  3 21:18:35 PST 2016
## Filename:      MazurkaPlugins/external/fftw/Makefile
## Syntax:        GNU Makefile
##
## Description:   This Makefile downloads and compiles the fft3 library.
##
## To run this makefile, type:
##     make everything
##

FFTW    = fftw-3.3.4
FFTWLIB = libfftw3.a
FFTWH   = fftw3.h
FFTWTAR = $(FFTW).tar.gz
FFTWURL = http://fftw.org/$(FFTWTAR)


all:
	@echo
	@echo "everything == do everything"
	@echo "download   == download fftw source code"
	@echo "unpack     == uncompress fftw source tarball"
	@echo "compile    == compile fftw library"
	@echo "install    == store compiled library in better place"
	@echo "clean      == delete unpacked source code directory"
	@echo "superclean == delete compiled library and source code"
	@echo


everything: download unpack compile install


download:
ifeq (,$(wildcard $(FFTWTAR)))
	wget $(FFTWURL)
endif


unpack:
ifeq (,$(wildcard $(FFTW)))
	tar xvzf $(FFTWTAR)
endif


compile:
ifeq (,$(wildcard lib/$(FFTWLIB)))
   ifeq (,$(wildcard $(FFTW)/Makefile))
	(cd $(FFTW); ./configure)
   endif
	(cd $(FFTW); make)
endif


install:
	mkdir -p include
	cp $(FFTW)/api/$(FFTWH) include/
	mkdir -p lib
	cp $(FFTW)/.libs/$(FFTWLIB) lib/


clean:
	-rm -rf $(FFTW)
	-rm $(FFTWTAR)


superclean: clean
	-rm -rf include
	-rm -rf lib



