####
Misc
####

* Optimisation

   Although the implementation likely has many optimisation opportunities, I respect and simply follow `this suggestion <https://en.wikiquote.org/wiki/Michael_A._Jackson#Principles_of_program_design,_1975>`_.
   Also I allow users to configure compiler optimisation options, and it may be possible to improve performance in certain cases by specifying link time optimisation or proper vectorisation.

* Absence of the final memory deallocations

   I intentionally omit the memory deallocation of the variables which survive throughout the simulation to simplify the whole things.
   On the other hand, I make sure that variables which are temporally used inside functions are properly freed to avoid memory leak.

* HTML theme of the documentation

   `Alabaster <https://alabaster.readthedocs.io/en/latest/>`_ is used after the colours are flipped.

