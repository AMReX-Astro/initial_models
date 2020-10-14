# nova

Options for init_1d are as follows:

## Gravity:

  * `g_type`: select gravity type; allowed values are:

    * `0`: constant gravity (calibrated for the base of the H1 layer)
    * `1`: 1/r^2 gravity
    * `2`: enclosed-mass gravity (calibrated for the
      base of H1 layer

  * `m`: mass of the star (in g)

## Core Profile:

  * `p_type`: select profile type; allowed values are:

    * `1`: isothermal profile
    * `2`: isentropic profile

  * `c_edge' select core profile edge; allowed values are:

    * `0`: extend initial model with profile
    * `1`: overwrite profile from base of envelope

## Smoothing:

  * `s_type`: select smoothing type; allowed values are:

    * `0`: no smoothing
    * `1`: Gaussian kernel smoothing
    * `2`: hyperbolic transition across the core-envelope boundary.
    * `3`: hyperbolic transition 0.5 slen below the core-envelope boundary
    * `4`: hyperbolic transition 2.0 slen below the core-envelope boundary
    * `5`: half-Gaussian transition below the core-envelope boundary
    * `6`: decaying-exponential transition below the core-envelope boundary

  * `s_length`: smoothing length
