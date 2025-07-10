PFiSM
=====

Time-evolution of simple phase-field model (PFM) coupled to temperature field.
PFM driving force is proportional to difference between local temperature and
melting temperature. Latent heat production from phase change affects temperature.
Alternatively, the temperature field can be fixed to a uniform value.
Boundary conditions can be periodic or Dirichlet or Neumann.

Implementation based on SAMRAI, Sundials and hypre libraries. Spatial discretization
is based on the finite volume approach. The time-integration is based on Sundials
and offers two options:
(i) a Backward-difference formula temporal discretization, using a Jacobian-Free
Newton Krylov (JNFK) solver based on Sundials CVODE functionalities to solve the
nonlinear problem in that implicit scheme. A multigrid preconditioner is implemented
using the hypre library to accelerate the solution of the linear systems within
the Newton iterations.
(ii) an IMEX (Implicit-Explicit) solver based on Sundials ARKODE functionalities.

Allows for testing and comparison of various time-integration schemes.

Conventions
-----------
Phase variable is 0 in liquid, 1 in solid.

Authors
-------

 * Jean-Luc Fattebert (fattebertj@ornl.gov)

Dependencies
------------

* [Hypre] (https://github.com/LLNL/hypre)

* [SAMRAI] (https://github.com/LLNL/SAMRAI)

* [Sundials] (https://github.com/LLNL/sundials)

* [HDF5] (https://support.hdfgroup.org/HDF5)
