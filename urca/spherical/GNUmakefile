NDEBUG := 
MPI    :=
OMP    :=

COMP := g++
FCOMP := gfortran

USE_REACT = TRUE
MKVERBOSE := t

SDC    :=

EBASE := initialmodel
# EOS and network
EOS_DIR := helmholtz

NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_reaclib/URCA-simple

Bpack := ./Make.package

EXTERN_SEARCH += . ../..

include $(MICROPHYSICS_HOME)/Make.Microphysics
