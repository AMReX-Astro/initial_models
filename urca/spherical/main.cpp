#include <new>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <time.h>

#include <microphysics.H>

#include <extern_parameters.H>

#include <fstream>

#include <network.H>
#include <eos.H>
#include <init_1d.H>

extern "C"
{
  void test_jacobian();
  void do_burn();
}


std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{

  //
  // Make sure to catch new failures.
  //
  amrex::Initialize(argc,argv);

  // save the inputs file name for later
  if (argc > 1) {
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }
  }


  // initialize the runtime parameters


  // we use a single file name for the extern name list and
  // the name list used by the initialization

  const int inputs_file_length = inputs_name.length();
  amrex::Vector<int> inputs_file_name(inputs_file_length);
  // initialize the runtime parameters

  for (int i = 0; i < inputs_file_length; i++) {
    inputs_file_name[i] = inputs_name[i];
  }

  runtime_init(inputs_file_name.dataPtr(), &inputs_file_length);

  init_extern_parameters();

  update_fortran_extern_after_cxx();


  // initialize Fortran Microphysics

  microphysics_initialize(small_temp, small_dens);


  // initialize C++ Microphysics

  eos_init(small_temp, small_dens);

  init_1d();

  amrex::Finalize();

  return 0;
}
