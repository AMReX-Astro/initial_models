[![github pages](https://github.com/AMReX-Astro/initial_models/actions/workflows/gh-pages.yml/badge.svg)](https://github.com/AMReX-Astro/initial_models/actions/workflows/gh-pages.yml) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8305500.svg)](https://doi.org/10.5281/zenodo.8305500)


# initial_models

These are various codes for generating 1-d HSE initial models for the
AMReX astrophysics codes.  We require amrex and StarKiller
Microphysics.

To build all of the dependencies, you need to set the `AMREX_HOME` and
`MICROPHYSICS_HOME` environment variables to point to their respective
top-level directories, e.g., as:

```
export AMREX_HOME=/path/to/amrex/
export MICROPHYSICS_HOME=/path/to/Microphysics/
```

where you replace `/path/to/` with the appropriate path on your system.

The routines are:

  * spherical

    generate an isentropic, self-gravitating WD model given a core
    temperature and density.


  * urca/spherical

    adaptation of the spherical initial model for a WD with the A=23
    Urca process. This initial model routine helps ensure the A=23
    Urca nuclei are in equilibrium in the initial spherically
    symmetric model.


  * test2

    generate an isentropic plane-parallel atmosphere with an entropy
    jump below to suppress convective overshoot.  This is used by the
    test2 and test_convect problems.


  * toy_atm

    similar to he_burn.  An isentropic layer is placed on top of an
    isothermal base.  A jump in temperature at the base of the
    isentropic layer is specified, with a linear transition between
    the base and isentropic layer.  The isentropic layer is continued
    down until reaching a cutoff temperature, at which point the model
    is isothermal.


  * lagrangian_planar

    This takes an existing initial model that is unequally gridded (in
    space) and maps it onto a uniform grid and puts it in HSE using
    our equation of state.  At the moment, it assumes a constant
    gravitational acceleration and plane-parallel.  The integration is
    down (up and down) from the location of the peak T in the model,
    and the temperature profile of the initial model is preserved.


  * kepler_hybrid

    These were used for our original set of wdconvect calculations.
    They take a 1-d model from the Kepler stellar evolution code
    (which may not be uniformly spaced), and put it into HSE on the
    MAESTRO grid.  This particular version forces the inner region of
    the star to be completely isentropic.
