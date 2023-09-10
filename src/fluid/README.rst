######
fluid/
######

Fluid solver.

* arrays

   Macros to access flow field.

* boundary_conditions

   Impose boundary conditions and exchange halo cells.

* compute_potential

   Poisson solver.

* compute_rhs

   Evaluate right-hand-side of the momentum equation.

* correct_velocity

   Project velocity from non-solenoidal to divergence-free field.

* decide_dt

   Decide time step size in the next step.

* init.c

   Allocator and flow-field loader.

* internal.h

   Private functions only used in this directory is declared.

* update_pressure.c

   Function to update pressure field.

* update_field

   Update velocity and temperature field using what is computed by ``compute_rhs``.

