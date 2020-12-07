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



extern "C"
{
  void test_jacobian();
  void do_burn();
  void do_initialization(const int inputs_name[], const int inputs_len);
  void init_1d();
  void init_1d_irreg();
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

  // we use a single file name for the extern name list and
  // the name list used by the initialization

  const int inputs_file_length = inputs_name.length();
  amrex::Vector<int> inputs_file_name(inputs_file_length);

  for (int i = 0; i < inputs_file_length; i++) {
    inputs_file_name[i] = inputs_name[i];
  }

  do_initialization(inputs_file_name.dataPtr(), inputs_file_length);

  init_1d();
  // init_1d_irreg();

  amrex::Finalize();

  return 0;
}
