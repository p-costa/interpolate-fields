#
# adapted from: https://fortran-lang.org/learn/building_programs/project_make
#
# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

.SUFFIXES:	.c .o .d .f90 .f .F90

# Project name
NAME := interpit

TARGET := $(NAME)

ROOTDIR := $(PWD)/..
SRCDIR := $(ROOTDIR)/src
APPDIR := $(ROOTDIR)/app
EXEDIR := $(ROOTDIR)/exe
INCDIR := $(SRCDIR)/include
CONFIGDIR := $(SRCDIR)/configs

EXE := $(EXEDIR)/$(TARGET)

# Configuration settings
FC := mpifort
FFLAGS :=
AR := ar rcs
LD := $(FC)
RM := rm -f
GD := ./.gen-deps.awk

include $(CONFIGDIR)/build.conf
include $(INCDIR)/flags.mk
include $(INCDIR)/compilers.mk

# List of all source files
SRCS := $(wildcard $(SRCDIR)/*.f90) $(wildcard $(APPDIR)/*.f90)
TEST_SRCS := 

# Define a map from each file name to its object file
obj = $(src).o
$(foreach src, $(SRCS) $(TEST_SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
DEPS := $(addsuffix .d, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
TEST_DEPS := $(addsuffix .d, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f90, %.exe, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
#all: $(LIB) $(TEST_EXE) $(EXE)
all: $(TEST_EXE) $(EXE)

# Create the static library from the object files
#$(LIB): #$(OBJS)
#	$(AR) $@ $^

# Link the test executables
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $^ -o $(EXE)
	mkdir -p $(EXEDIR)/data

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: % | %.d
	$(FC) $(FFLAGS) -c -o $@ $<

# Process the Fortran source for module dependencies
$(DEPS) $(TEST_DEPS): %.d: %
	$(GD) $< > $@

# Define all module interdependencies
include $(DEPS) $(TEST_DEPS)
$(foreach dep, $(OBJS) $(TEST_OBJS), $(eval $(dep): $($(dep))))

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(filter %.d, $(DEPS) $(TEST_DEPS)) $(filter %.exe, $(TEST_EXE)) $(wildcard *.mod) $(EXE)
