PFiSM
=====

Implementation of simple phase-field model (PFM) based on SAMRAI, Sundials
and hypre. Includes coupled solver for temperature field.
PFM driving force is proportional to difference between local temperature
and melting temperature. Latent heat production from phase change affects
temperature.

Authors
-------

 * Jean-Luc Fattebert (fattebertj@ornl.gov)

Dependencies
------------

* [Hypre] (https://github.com/LLNL/hypre)

* [SAMRAI] (https://github.com/LLNL/SAMRAI)

* [Sundials] (https://github.com/LLNL/sundials)

* [HDF5] (https://support.hdfgroup.org/HDF5)
