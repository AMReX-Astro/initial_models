NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t


# EOS and network
EOS_DIR := gamma_law_general

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/H_He.net


# the local sources needed to build the initial model generator
MODEL_SOURCES := init_1d.f90


include ../GInitModels.mak


