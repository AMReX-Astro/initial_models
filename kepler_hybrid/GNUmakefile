# Makefile for the init_1d initial model generator

NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# EOS and network
EOS_DIR := helmholtz

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/ignition_wdconvect.net


# the local sources needed to build the initial model generator
# MODEL_SOURCES := init_1d.f90
MODEL_SOURCES := init_1d_irreg.f90


include ../GInitModels.mak


