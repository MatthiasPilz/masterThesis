###############################################################################
# Set different path variables                                                #
###############################################################################

# Working directory
TOP := $(shell pwd)

# Executable
EXE    = run_ADER
EXEDIR = $(TOP)/../exe

#Main directory
MAIN   = main

#Misc directory
MSC   =  $(TOP)/../msc

INCDIR = $(TOP)/main
LIBDIR = $(TOP)/../lib
OBJDIR = $(LIBDIR)/obj
MATHDIR= $(TOP)/../math/
CONVDIR= $(TOP)/../conv/


###############################################################################
# Set compiler and compiler flags                                             #
###############################################################################

# Compiler name e.g. gcc or icc
CC = gcc

# BLAS resp GSL
CFLAGS += -g
libsys += -lgsl -lgslcblas

# standard math library, sometimes has to go last
libsys += -lm

# optimization, e.g. -g, -O, -O2, -O3
OFLAGS = -O2
OFLAGS += -g

# warnings, e.g. -Wall
WARN   =  #-Wall

CFLAGS += $(OFLAGS) -I$(INCDIR) $(WARN)
CFLAGS := $(strip $(CFLAGS))


