# A set of useful macros for putting together one of the initial model
# generator routines

# include the main Makefile stuff
include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# default target (make just takes the one that appears first)
ALL: init_1d.$(suf).exe


#-----------------------------------------------------------------------------
# core FBoxLib directories
FBOXLIB_CORE := Src/BaseLib
FPP_DEFINES += -DAMREX_DEVICE=""


#-----------------------------------------------------------------------------
# Microphysics directories needed

EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks

# the helmeos has an include file -- also add a target to link the table
# into the problem directory.
ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
  ALL: table
  EOS_PATH := $(EOS_TOP_DIR)/helmholtz
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi


MICROPHYS_CORE := $(MICROPHYSICS_HOME)/EOS \
                  $(MICROPHYSICS_HOME)/networks \
                  $(MICROPHYSICS_HOME)/interfaces \
                  $(MICROPHYSICS_HOME)/unit_test


# add in the network, EOS
MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                  $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \

# get any additional network dependencies
include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES


# explicitly add in any source defined in the build directory
f90sources += $(MODEL_SOURCES)
f90sources += $(READER_SOURCES)


#-----------------------------------------------------------------------------
# core FBoxLib directories
Fmpack := $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir))
Fmincs :=

# auxillary directories
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)



# we need a probin.f90, since the various microphysics routines can
# have runtime parameters
F90sources += probin.F90

PROBIN_TEMPLATE := $(MICROPHYSICS_HOME)/unit_test/dummy.probin.template
PROBIN_PARAMETER_DIRS = $(INITIAL_MODEL_HOME)
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE) $(NETWORK_TOP_DIR)


PROBIN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.F90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.F90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "

read_mesa.f90: $(INITIAL_MODEL_HOME)/read_mesa/read_mesa.f90
	cp ../read_mesa/read_mesa.f90 .
	# $(LINK.f90) -o read_mesa.o $(objects) $(libraries)
	@echo SUCCESS






# vpath defines the directories to search for the source files

#  VPATH_LOCATIONS to first search in the problem directory
#  Note: GMakerules.mak will include '.' at the start of the
VPATH_LOCATIONS += . $(Fmlocs)



# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90:
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/makebuildinfo.py \
           --modules "$(Fmdirs) $(MICROPHYS_CORE) $(UNIT_DIR)" \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --fboxlib_home "$(FBOXLIB_HOME)" \
           --source_home "$(MICROPHYSICS_HOME)" \
           --network "$(NETWORK_DIR)" \
           --integrator "$(INTEGRATOR_DIR)" \
           --eos "$(EOS_DIR)"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90

init_1d.$(suf).exe: $(objects)
	$(LINK.f90) -o init_1d.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

# include the fParallel Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean::
	$(RM) probin.F90
	$(RM) build_info.f90
	$(RM) read_mesa.f90
