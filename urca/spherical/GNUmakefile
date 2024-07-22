PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

COMP := g++
FCOMP := gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

BL_NO_FORT = TRUE

MAX_ZONES := 4096

DEFINES += -DNPTS_MODEL=$(MAX_ZONES)
SDC    :=

EBASE := initialmodel
# EOS and network
EOS_DIR := helmholtz

NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_reaclib/URCA-simple

Bpack := ./Make.package
Blocs := . ../..

EXTERN_SEARCH += . ../..

include $(MICROPHYSICS_HOME)/Make.Microphysics
