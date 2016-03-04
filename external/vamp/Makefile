##
## Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
## Creation Date: Thu Mar  3 21:31:08 PST 2016
## Last Modified: Thu Mar  3 21:31:11 PST 2016
## Filename:      MazurkaPlugins/external/vamp/Makefile
## Syntax:        GNU Makefile
##
## Description:   This Makefile downloads and compiles the vamp sdk.
##
## To run this makefile, type:
##     make everything
##

VAMP    = vamp-plugin-sdk-2.6
VAMPTAR = $(VAMP).tar.gz
VAMPURL = https://code.soundsoftware.ac.uk/attachments/download/1520/$(VAMPTAR)
VAMPLIB = libvamp-sdk.a


all:
	@echo
	@echo "everything == do everything"
	@echo "download   == download vamp source code"
	@echo "unpack     == uncompress fftw source tarball"
	@echo "compile    == compile fftw library"
	@echo "install    == setup interface directories"
	@echo "clean      == delete unpacked source code directory"
	@echo "superclean == delete compiled library and source code"
	@echo


everything: download unpack compile install fix


download:
ifeq (,$(wildcard $(VAMPTAR)))
	wget $(VAMPURL)
endif


unpack:
ifeq (,$(wildcard $(VAMP)))
	tar xvzf $(VAMPTAR)
endif


compile:
ifeq (,$(wildcard lib/$(VAMPLIB)))
	(cd $(VAMP); ./configure)
   ifeq (,$(wildcard $(VAMP)/Makefile))
	# The compile will fail due to misisng libsndlib.a
	# but that is OK, since only interested in static library
	# for Vamp SDK.
	-(cd $(VAMP); make -f build/Makefile.osx)
   endif
endif


install:
	mkdir -p include
	cp $(VAMP)/vamp/vamp.h include/
	cp $(VAMP)/vamp-sdk/*.h include/
	mkdir -p lib
	cp $(VAMP)/$(VAMPLIB) lib/


fix:
	sed 's/<vamp\/vamp.h>/"vamp.h"/' include/PluginAdapter.h > temp.h
	mv -f temp.h include/PluginAdapter.h


clean:
	-rm -rf $(VAMP)
	-rm $(VAMPTAR)


superclean: clean
	-rm -rf include
	-rm -rf lib



