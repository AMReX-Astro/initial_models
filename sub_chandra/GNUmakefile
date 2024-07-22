PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 1

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

BL_NO_FORT = TRUE

MAX_ZONES := 16384

DEFINES += -DNPTS_MODEL=$(MAX_ZONES)

# programs to be compiled
EBASE := initialmodel

# EOS and network
EOS_DIR := helmholtz

NETWORK_DIR := subch_simple

Bpack := ./Make.package
Blocs := . ..

EXTERN_SEARCH += . ..

include $(MICROPHYSICS_HOME)/Make.Microphysics
