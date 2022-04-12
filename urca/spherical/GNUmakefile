PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

COMP := g++
FCOMP := gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE
USE_CXX_EOS = TRUE

MAX_ZONES := 16384

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
