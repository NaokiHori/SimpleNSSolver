####
src/
####

All functions are implemented here.
For the directories, see the corresponding README file.

* array.c

   Function to initialise, destruct, load, and save multi-dimensional arrays including halo cells are implemented.

* config.c

   Environment variable loader which is called when the solver is launched to acquire the runtime-parameters by the user is implemented.

* domain.c

   Information about the coordinate systems, the grid configuration, and the domain parallelisation is included.

* linear_system.c

   Functions to initialise and destruct the tri-diagonal linear system solver are included.

* main.c

   Main function is here.

* memory.c

   Utility functions and global parameters are defined.

* runge_kutta.c

   Runge-Kutta coefficients are defined.

* save.c

   Functions to save flow field and parameters for restart are included.

* statistics.c

   Routines to collect statistical data are implemented.

* tdm.c

   Kernel functions to solve tri-diagonal matrices are implemented.

* timer.c

   A function to obtain the current wall time is implemented.

