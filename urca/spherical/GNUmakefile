NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

SDC    :=

# EOS and network
EOS_DIR := helmholtz

NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_reaclib/URCA-simple

# the local sources needed to build the initial model generator
MODEL_SOURCES := init_1d.f90
MODEL_SOURCES += urca_composition.f90

include GInitModels.mak
