************
Introduction
************

This repository holds the routines that are used to make hydrostatic
initial models for the AMReX Astrophysics Suite codes MAESTROeX and
Castro.  It is designed to work with the equations of state and
reaction networks defined in the StarKiller Microphysics repository.

There are several routines that are roughly broken into two categories:
those that create simple paramaterized models on their own and those
that adjust models from stellar evolution codes.

All of these model builders enforce a discrete form of hydrostatic equilibrium:

.. math::

   \frac{p_{i+1} - p_i}{\Delta r} = \frac{1}{2} (\rho_i + \rho_{i+1} ) g_{i+1/2}

Together with the equation of state.  We typically integrate outward from the center
of the star or base of an atmosphere, using the known :math:`i` state to find the
:math:`i+1` state above:

.. math::

   p_{i+1} = p_i + \frac{\Delta r}{2} (\rho_i + \rho_{i+1} ) g_{i+1/2}

One additional equation is
needed, either a specified temperature profile (or isothermal) or
isentropic, or some mix of these in different regions in the star.

If we are constraining to a temperature :math:`T_c`, then we have an implicit
equation for :math:`\rho_{i+1}` to solve:

.. math::

   p(\rho_{i+1}, T_c) = p_i + \frac{\Delta r}{2} (\rho_i + \rho_{i+1} ) g_{i+1/2}

And similar for a constant entropy.  

Simple parameterized models
---------------------------

  * ``low_mass_convective_star``

    Generate an spherical, fully convective star of specified mass.
    This is meant to be used with the ``gamma_law`` equation of state
    and essentially generates a polytrope.

  * ``spherical``

    Generate an isentropic, self-gravitating WD model given a core
    temperature and density.  This is generally used for making a
    white dwarf.


  * ``sub_chandra``

    Create a sub-Chandra C/O white dwarf with a thin He layer on the surface.  Both the
    mass of the WD and the He layer are inputs parameters.

  * ``test2``

    Generate an isentropic plane-parallel atmosphere with an entropy
    jump below to suppress convective overshoot.  This is used by the
    ``reacting_bubble`` (originally called "test2") and
    ``test_convect`` problems.


  * ``toy_atm``

    An isentropic layer is placed on top of an isothermal base.  A
    jump in temperature at the base of the isentropic layer is
    specified, with a linear transition between the base and
    isentropic layer.  The isentropic layer is continued down until
    reaching a cutoff temperature, at which point the model is
    isothermal.  This was used for some of the X-ray burst models.


stellar evolution code converters
---------------------------------

  * ``lagrangian_planar``

    This takes an existing initial model that is unequally gridded (in
    space) and maps it onto a uniform grid and puts it in HSE using
    our equation of state.  At the moment, it assumes a constant
    gravitational acceleration and plane-parallel.  The integration is
    down (up and down) from the location of the peak T in the model,
    and the temperature profile of the initial model is preserved.


  * ``kepler_hybrid``

    These were used for our original set of wdconvect calculations.
    They take a 1-d model from the Kepler stellar evolution code
    (which may not be uniformly spaced), and put it into HSE on the
    MAESTRO grid.  This particular version forces the inner region of
    the star to be completely isentropic.

  * ``massive_star``

    This takes a model for a massive star and puts into to HSE.  It
    uses the temperature structure from the model and Ye or X,
    depending on whether the state is in nuclear statistical
    equilibrium.
