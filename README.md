# PHENOM 

[![](https://github.com/twguest/phenom/actions/workflows/testing.yml/badge.svg)]( https://github.com/twguest/phenom/actions/workflows/testing.yml)[![](https://img.shields.io/pypi/v/phenom.svg)](https://pypi.python.org/pypi/phenom_xfel)


A **phenom**enological model of X-ray Free Electron Laser (XFEL) radiation.

The PHENOM python package is designed to provide a simple, robust and computationally efficient method for generating representations of the complex wavefield of X-ray Free Electron Laser pulses. By making use of approximate representations of pulse wavefront and [spectra](https://www.osapublishing.org/abstract.cfm?URI=ol-35-20-3441), phenom allows large ensembles of photon pulses with arbitrary statistics to be generated in a truly python-ised manner.

## Getting Started
At the command line::

    $ pip install phenom-xfel

To check that your instillation has worked, open iPython and try::

    $ import phenom
    
For more information on usage, see the docs [here](https://twguest.github.io/phenom).

## Features
The phenomenological SASE model consists of two primary functions:

1. SASE pulse spectrum
2. Gaussian beam profile

Application of the model is based on the assumption that time-varying fluctuations in these wavefield components are sufficient to describe much of the shot-to-shot statistical properties of XFEL radiation.

## Examples
Phenom has been designed to require minimal knowledge of the XFEL process prior to generating your first pulse.

Details on generating these pulses can be found in the [documentation](https://twguest.github.io/phenom).

We highlight three

## Special Thanks:

[1]: https://github.com/JunCEEE
[![github]([https://cloud.githubusercontent.com/assets/17016297/18839843/0e06a67a-83d2-11e6-993a-b35a182500e0.png](https://avatars.githubusercontent.com/u/3163505?v=4)https://avatars.githubusercontent.com/u/3163505?v=4)][1]
