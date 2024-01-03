ESCAPE_HOME := .

OBJECTS := Graph.o GraphIO.o

TARGETS := libescape.a

all: clustering

libescape.a : $(OBJECTS)
	ar cruv $@ $^

include common.mk

clustering: libescape.a
	$(MAKE) -C clustering

cleanclustering:
	$(MAKE) -C clustering clean cleandep

.PHONY: clustering



