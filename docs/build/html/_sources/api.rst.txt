===
API
===

SASE Source
###########

The simplest way to generate represenations of the SASE field is using the SASE
wavefront function. 

The properties of this wavefield are based on a model of XFEL radiation containing 3 components:

A complex wavefront
=======================

.. autofunction:: phenom.gaussian.complex_gaussian_beam

A complex time spectrum
============================

.. autofunction:: phenom.spectrum.linear_SASE_spectrum

A complex pointing angle tilt
==================================
.. autofunction:: phenom.wavefront_tools.wavefront_tilt
 
These are combined in the sase pules module.

.. autofunction:: phenom.source.sase_pulse

Generating a SASE pulse looks something like this:


.. code-block:: python

    import numpy as np
    from phenom.source import sase_pulse
    
    x = np.linspace(-500e-06, 500e-06, 512)
    y = np.linspace(-500e-063, 500e-06, 512)
    t = np.linspace(-100e-15, 100e-15, 512)
    
    electric_field  = sase_pulse(x = x, 
                                  y = y,
                                  t = t,
                                  photon_energy = 9200.,
                                  pulse_energy = 1e-03,
                                  pulse_duration = 25e-15,
                                  bandwidth = 1e03,
                                  sigma = 100e-06,
                                  div = 1e-03,
                                  x0 = 0.,
                                  y0 = 0.,
                                  t0 = 2e-15,
                                  theta_x = 0,
                                  theta_y = 0.                                                              