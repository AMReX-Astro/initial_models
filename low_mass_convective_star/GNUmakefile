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
EOS_DIR := gamma_law_general

NETWORK_DIR := general_null
NETWORK_INPUTS := H_He.net

Bpack := ./Make.package
Blocs := .

EXTERN_SEARCH += . ..

include $(MICROPHYSICS_HOME)/Make.Microphysics
