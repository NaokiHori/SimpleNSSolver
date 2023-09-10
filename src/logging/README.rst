########
logging/
########

Functions to monitor the running simulation.

* divergence.c

   Compute and output maximum local divergence of the flow field.

* energy.c

   Compute and output the total squared velocity (in each direction) and the total squared temperature of the flow field.

* internal.h

   Private functions which are only used by this directory is declared.

* main.c

   Call other logging functions.

* momentum.c

   Compute and output the net momentum in each direction.

