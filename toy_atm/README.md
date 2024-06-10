# `toy_atm`

Create a 1-d hydrostatic, atmosphere with an isothermal region
(`T_star`) representing the NS, a hyperbolic tangent rise to a
peak temperature (`T_base`) representing the base of an accreted
layer, an isoentropic profile down to a lower temperature (`T_lo`),
and then isothermal. This can serve as an initial model for a
nova or XRB.

The temperature profile is:

```
        ^
        |
 T_base +        /\
        |       /  \
        |      /  . \
 T_star +-----+      \
        |     .   .   \
        |              \
        |     .   .     \
 T_lo   +                +-----------
        |     .   .
        +-----+---+---------------> r
        |      \  /
        |      delta
        |< H_star>|
```

We take `dens_base`, the density at the base of the isentropic layer
as input.  The composition is "ash" in the lower isothermal region
and "fuel" in the isentropic and upper isothermal regions.  In the
linear transition region, we linearly interpolate the composition.

The fuel and ash compositions are specified by the `fuel?_name`,
`fuel?_frac` and `ash?_name`, `ash?_frac` parameters (name of the species
and mass fraction).

The model is placed into HSE by the following differencing:

```
   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
```

This will be iterated over in tandem with the EOS call,

```
   P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
```

## Parameters

The following parameters can be set in the inputs file (prefixed with
`problem.`).

### Thermodynamics

* `dens_base` : density at the base of the fuel layer

* `T_star` : temperature in the isothermal region beneath the fuel layer

* `T_base` : temperature at the base of the fuel layer

* `T_lo` : lowest temperature in the fuel layer.  Below this we switch from
  being isentropic to isothermal

* `H_star` : Height to base of the fuel layer

* `delta` : with the transition region between `T_star` and `T_base`

* `low_density_cutoff` : the density below which we don't worry about HSE.

### Composition:

The composition of the fuel is set by a pair of parameters, one giving
the name of the species and the other giving the mass fraction.  These
all have the form:

* `fuel?_name`, `fuel?_frac`, where `?` can be 1, 2, ... 7.  At least the
  first needs to be set.

The ash composition is similarly set via:

* `ash?_name`, `ash?_frac`, where `?` can be 1, 2, ... 7.  At least the
  first needs to be set.

Finally:

* `smallx` : the smallest allowed mass fraction

### Domain

The domain size and resolution are set via:

* `xmin` : coordinate of bottom of the model

* `xmax` : coordinate of the top of the model

* `nx` : number of points.

These should be chosen to ensure that the grid spacing, dx = (xmax - xmin) / nx
corresponds to the finest resolution in the simulation you will be running.

### Gravity

There are 2 options for gravity: constant and a `1/r**2` profile:

* `g_const` : the constant gravitational acceleration

* `M_enclosed` : for `1/r**2` gravity, the mass enclosed beneath the
  model.  This will be used together with `xmin` to compute the
  gravitational acceleration at the base of the model.

* `do_invsq_grav` : do we use constant? or `1/r**2` gravity?


## Science problems where this is used

This setup is used for ``toy_convect`` and for the many of the XRB problems.
The init code takes an input file:

* `inputs.xrb_mixed.hi_dens.tall.CNO`

   This is the `_params` file that is used for the `xrb_mixed` file

* `inputs.xrb_mixed.hi_dens.tall2.CNO`

   This version differs from the above one with more buffer beneath
   the convective layer


