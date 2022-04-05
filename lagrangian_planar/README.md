# lagrangian_planar

This is a generic initial model routine that takes an starting model
from a 1-d plane-parallel Lagrangian stellar evolution code and puts
it onto a uniform grid and then enforces HSE.  This is designed for
an X-ray burst (for example).

The gravity can be constant or a 1/r**2 gravity.

The HSE integration is done from the index of the peak temperature,
integrating up and down from there, using the interpolated temperature
from the initial model and a discrete form of HSE.
