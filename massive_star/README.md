# massive star initial model

This is the initial model routine for the massive star problem.

The file: `15m_500_sec.txt` is a 15 solar mass MESA initial model, with
an aprox21 composition.

It can be converted into the subset of nuclei for aprox19 via
`convert_21_to_19.py`

This setup uses Ye as the primary composition variable from the initial
model in regions that are in NSE.

Note: you should ensure that the NSE conditions in the inputs file match
those of your simulation, so the model will be properly in HSE.

Also note that when running with > 16384 zones, you need to do:

```
ulimit -s 32768
```

since the arrays are put on the stack.

## Overall algorithm

The basic HSE algorithm proceeds as:

* Map the raw MESA model onto a uniform grid

* Reset the composition of any zones that are in NSE using the NSE
  table

* (optional) If the first MESA r is at a larger radius than the
  innermost uniform model zone, integrate inward from the first MESA
  point at constant entropy to find the central density in HSE at our
  model resolution.  This is only necessary if the MESA model is very
  low resolution.

* Integrate outward from the center taking T and composition from the
  MESA model and enforcing HSE and NSE with our EOS.
