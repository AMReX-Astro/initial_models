PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 1

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

# programs to be compiled
EBASE := initialmodel

# EOS and network
EOS_DIR := helmholtz

NETWORK_DIR := aprox19
USE_NSE := TRUE

Bpack := ./Make.package ../read_mesa/Make.package
Blocs := . ../read_mesa

EXTERN_SEARCH += . ..

include $(MICROPHYSICS_HOME)/Make.Microphysics
