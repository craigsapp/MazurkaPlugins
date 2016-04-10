

BASE   := $(shell basename `pwd`)
BINDIR = ../bin

all: index htmlcode


index:
	../bin/makepluginindex ../template.html ./


htmlcode:
	
	../bin/codecommentary src/$(BASE).cpp > src/$(BASE).cpp.html
	../bin/codecommentary src/$(BASE).h > src/$(BASE).h.html




