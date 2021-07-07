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

NETWORK_DIR := subch
#NETWORK_INPUTS := triple_alpha_plus_o.net

Bpack := ./Make.package
Blocs := .

EXTERN_SEARCH += . ..

include $(MICROPHYSICS_HOME)/Make.Microphysics
