ESCAPE_HOME := ../

TARGETS := RTRex nucleus

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

count_% : count_%.o ../libescape.a
	$(COMPILE_AND_LINK)

