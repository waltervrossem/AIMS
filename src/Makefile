#
# Author: Daniel R. Reese <daniel.reese@obspm.fr>
# Copyright (C) Daniel R. Reese and contributors
#
# 2020/05/08 Added type checking, linting, and unit test framework.
#            Steven J. Hale <s.j.hale@bham.ac.uk>
#
################################################################################

ROOT := ..
BINDIR := $(ROOT)/bin
ETCDIR := $(ROOT)/etc
SRCDIR := $(ROOT)/src

F90 := gfortran
FLG := -O2 #-fcheck=all #-fbounds-check
FPY := f2py
#FPY := f2py3 # for python3
#FPY := python3 /full/address/to/f2py # another way of working with python3

################################################################################

.PHONY: all build lint test doccov clean

all: build

################################################################################

##
# Tell the user what the targets are.
##

help:
	@echo
	@echo "Targets:"
#	@echo "   build      Type-check and compile all modules."
	@echo "   build      Compile all modules."
	@echo "   lint       Check that all modules satisfy the PEP 8 coding standard."
	@echo "   test       Run all unit tests."
	@echo "   doccov     Calculate documentation coverage."
	@echo "   clean      Remove unversioned files."
	@echo

################################################################################

# Run in --silent mode unless the user sets VERBOSE=1 on the
# command-line.

ifndef VERBOSE
.SILENT:
endif

ifdef REPORTS
PYLINTFLAGS = --reports=y
endif

################################################################################

type:
	PYTHONPATH=$(SRCDIR) \
	MYPYPATH=$(ETCDIR)/Mypyrc/stubs \
	python3 -m mypy $(SRCDIR) --config-file=$(ETCDIR)/Mypyrc/mypy.ini

compile:
	$(FPY) -c aims_fortran.f90 -m aims_fortran --fcompiler=$(F90) --f90flags=$(FLG)
	find $(SRCDIR) -name "*.py"  -exec python3 -m py_compile '{}' \;

#build: type compile
build: compile

lint:
	PYTHONPATH=$(SRCDIR) \
	PYLINTHOME=$(ETCDIR)/Pylintrc \
	PYLINTRC=$(ETCDIR)/Pylintrc/pylintrc \
	python3 -m pylint $(PYLINTFLAGS) $(SRCDIR)/*.py

test:
	PYTHONPATH=$(SRCDIR) \
	python3 -m pytest --cov-config=$(ETCDIR)/Coveragerc/coveragerc --cov -v -ra -s

doccov:
	PYTHONPATH=$(SRCDIR) \
	docstr-coverage $(SRCDIR)

clean :
	git clean -dxf

################################################################################
