#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = FindNeutrinos_FiltersDict
SHLIB = libFindNeutrinos_Filters.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell basictool-config --includes)
INCFLAGS += -I$(OPENCV_INCDIR)

LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell basictool-config --libs)
LDFLAGS += -L$(OPENCV_LIBDIR) -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE
