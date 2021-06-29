******************
Runtime Parameters
******************

Microphysics/extern parameter format
------------------------------------

The microphysics/extern parameter definitions take the form of:

::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the priority is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.

The documentation below for the Castro C parameters is
automatically generated, using the comments in the \_cpp_parameters
file.


.. toctree::

   runtime_parameters
