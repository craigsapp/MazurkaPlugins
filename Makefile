##
## Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
## Creation Date: Fri Mar  4 04:15:22 PST 2016
## Last Modified: Fri Mar  4 04:15:26 PST 2016
## Filename:      MazurkaPlugins/external/Makefile
## Syntax:        GNU Makefile
##
## Description:   This Makefile prepares the external libraries and
##                compiles the target plugins
##
## To run this makefile, type:
##     make everything
##

.PHONY: external

all: plugin


plugin:
	$(MAKE) -f Makefile.plugin


install:
	mkdir -p ~/Library/Audio/Plug-Ins/Vamp
	cp lib-osx/mazurka-plugins.dylib ~/Library/Audio/Plug-Ins/Vamp/


external:
	(cd external; $(MAKE) everything)


clean:		
	-rm -rf obj-osx
	(cd external; $(MAKE) clean)


superclean: clean
	-rm -rf lib-osx
	(cd external; $(MAKE) superclean)










