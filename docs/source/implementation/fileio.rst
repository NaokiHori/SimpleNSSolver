
.. _fileio:

###############
`src/fileio.c`_
###############

.. _src/fileio.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fileio.c

``fileio.c`` contains functions that handle file input/output operations.
Mostly they are designed to read and write files in `NPY format <https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html>`_.

The reasons for using the ``NPY`` format in this project are as follows:

#. Being binary

   I should avoid saving arrays in non-binary (text, e.g. ASCII) format, as the original floating-point data will be modulated and the file size increases even further.

#. Being simple

   To be able to recover the original data later, I need to store the datatype and shape as metadata with the raw data.
   In a ``NPY`` file, this information is stored in the header as ``descr`` and ``shape`` values respectively.

   .. note::

      The `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_ offers the same or more general functionality.
      However, it is too general to store simple scalars, vectors, and structured arrays used in this library (see `this page <https://numpy.org/doc/1.13/neps/npy-format.html#alternatives>`_).
      Moreover, I would like to reduce the number of external libraries to keep things simple.

#. Being general

   `NumPy <https://numpy.org>`_ is one of the most popular libraries for scientific computing (especially for post-processing and visualisation with `Matplotlib <https://matplotlib.org/stable/>`_), and the ``NPY`` format is its native binary format.
   Thanks to its simple structure, this format is easily accessible from other languages (e.g. `Rust <https://docs.rs/npy/latest/npy/>`_, `MATLAB <https://github.com/kwikteam/npy-matlab>`_).

.. seealso::

   Reading and writing ``NPY`` headers rely on my library `SimpleNpyIO <https://github.com/NaokiHori/SimpleNpyIO>`_.

