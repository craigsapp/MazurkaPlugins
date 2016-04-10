
BIN = bin2

.PHONY	: download

all: pages index download


download:
	(cd download; make)

index:
	$(BIN)/makecommandindex index.con > index.html

pages:
	$(BIN)/makepages .






