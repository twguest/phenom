# phenom

[![](https://github.com/twguest/phenom/actions/workflows/testing.yml/badge.svg)]( https://github.com/twguest/phenom/actions/workflows/testing.yml)[![](https://img.shields.io/pypi/v/phenom.svg)](https://pypi.python.org/pypi/phenom_xfel)


A **phenom**enological model of X-ray Free Electron Laser (XFEL) radiation.

The PHENOM python package is designed to provide a simple, robust and computationally efficient method for generating representations of the complex wavefield of X-ray Free Electron Laser pulses. By making use of approximate representations of pulse wavefront and [spectra](https://www.osapublishing.org/abstract.cfm?URI=ol-35-20-3441), phenom allows large ensembles of photon pulses with arbitrary statistics to be generated in a truly python-ised manner.

## Getting Started
At the command line::

    $ pip install phenom-xfel

To check that your instillation has worked, open iPython and try::

    $ import phenom
    
## Examples
Phenom has been designed to require minimal knowledge of the XFEL process prior to generating your first pulse.

1. [Getting Started](https://twguest.github.io/phenom/notebooks/sase_model_pt1.html)
2. [Tutorials](https://twguest.github.io/phenom/notebooks/sase_model_pt2.html).
3. [Integrating with WPG](https://twguest.github.io/phenom/notebooks/phenom_to_wpg.html).

More details on generating these pulses can be found in the [documentation](https://twguest.github.io/phenom).

