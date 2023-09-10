#######################
Nusselt number relation
#######################

The Nusselt number in the general coordinate system is defined as

.. math::

   Nu
   \equiv
   - \ave{\vat{\der{T}{x}}{\text{wall}}}{y,t}
   =
   - \ave{\frac{\int \vat{\der{T}{x}}{\text{wall}} dy}{\int dy}}{t}
   =
   - \ave{\frac{\int \der{\gx}{x} \vat{\der{T}{\gx}}{\text{wall}} \der{y}{\gy} d\gy}{\int \der{y}{\gy} d\gy}}{t}.

After discretised, I have

.. math::

   - \ave{\frac{\sum_j \frac{1}{\Delta x_{\text{wall}}} \vat{\diffe{T}{\gx}}{\text{wall}} \Delta y_j }{\sum_j \Delta y_j}}{t}
   =
   - \frac{1}{l_y} \ave{\sum_j \frac{1}{\Delta x_{\text{wall}}} \vat{\diffe{T}{\gx}}{\text{wall}} \Delta y_j}{t},

where :math:`l_y` is the domain size in the :math:`y` direction.

In this part, I look into the details of the additional relations which are satisfied in the Rayleigh-Bénard convections in this project.

*********
Heat flux
*********

.. include:: heat_flux.rst

****************
Energy injection
****************

.. include:: injection.rst

**************************
Kinetic energy dissipation
**************************

.. include:: kinetic.rst

**************************
Thermal energy dissipation
**************************

.. include:: thermal.rst

