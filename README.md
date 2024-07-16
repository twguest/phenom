# phenom

![pre-commit](https://github.com/twguest/phenom/actions/workflows/pre-commit.yml/badge.svg)
![docs](https://github.com/twguest/phenom/actions/workflows/static.yml/badge.svg)
[![](https://img.shields.io/pypi/v/phenom.svg)](https://pypi.python.org/pypi/phenom_xfel)
[![Documentation](https://img.shields.io/badge/documentation-online-blue)](https://twguest.github.io/phenom/)
![Language](https://img.shields.io/badge/language-python-blue)

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

## Citation:
The use of this package and the methods applied therein should be acknowledged using the following citation:
- Guest, T. W., R. Bean, R. Kammering, G. van Riessen, A. P. Mancuso, and B. Abbey. “A Phenomenological Model of the X-Ray Pulse Statistics of a High-Repetition-Rate X-Ray Free-Electron Laser.” IUCrJ 10, no. 6 (November 1, 2023). https://doi.org/10.1107/S2052252523008242.

## Cited By:
Below is a list of research which have acknowledged the application of this model
- E, Juncheng, Carsten Fortmann-Grote, Trey Guest, Egor Sobolev, Luca Gelisio, Richard Bean, and Adrian P. Mancuso. “SimEx-Lite: Easy Access to Start-to-End Simulation for Experiments at Advanced Light Sources.” In Advances in Computational Methods for X-Ray Optics VI, edited by Oleg Chubar and Takashi Tanaka, 22. San Diego, United States: SPIE, 2023. https://doi.org/10.1117/12.2677299.


