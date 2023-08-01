#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:01:52 2023

@author: twguest
"""

import h5py as h5
import numpy as np

from phenom.utils import e2wav


def wavefront_tilt(x, y, theta_x, theta_y, photon_energy):
    """
    Define a two-dimensional complex tilted plane-wave, tilted along transverse pointing vectors kx and ky.

    :param x: numpy.ndarray
        A 1D array representing the x-coordinates. The array size corresponds to nx.

    :param y: numpy.ndarray
        A 1D array representing the y-coordinates. The array size corresponds to ny.

    :param theta_x: float
        Horizontal tilting angle in radians, defining the horizontal pointing vector. Range: [-pi, pi].

    :param theta_y: float
        Vertical tilting angle in radians, defining the vertical pointing vector. Range: [-pi, pi].

    :param photon_energy: float
        The photon energy to be used in the wavefront calculation.

    :return: numpy.ndarray
        A 2D complex array representing the tilted plane-wave.

    :Note:
        The mesh is generated from the 'x' and 'y' parameters (using numpy's meshgrid), and kx and ky are computed
        based on the tilting angles 'theta_x' and 'theta_y' and the photon energy.
        The 'e2wav' function is assumed to be used for energy to wavelength conversion within the function.
    """
    xx, yy = np.meshgrid(x, y)
    kx = 2 * np.pi * np.sin(theta_x) / e2wav(photon_energy)
    ky = 2 * np.pi * np.sin(theta_y) / e2wav(photon_energy)

    tilt = np.exp(1j * (kx * xx + ky * yy))
    tilt /= np.sqrt(np.sum(abs(tilt**2)))
    return tilt


def geometric_focus(fwhm, divergence):
    """
    Calculate the focal length of a converging beam of size fwhm, approximating that we are viewing the beam far
    downstream of the beam waist.

    :param fwhm: float
        The full width at half maximum (FWHM) of the converging beam.

    :param divergence: float
        The divergence angle of the beam in radians, typically a small angle.

    :return: float
        The focal length of the converging beam.

    :Note:
        This function uses a simple geometrical approximation considering the beam size and divergence,
        and is most suitable for scenarios where the observation point is far from the beam waist.
    """
    return fwhm / np.tan(divergence)
