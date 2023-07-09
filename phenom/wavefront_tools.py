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
    define a two-dimensional complex tilted plane-wave, tilted along transverse pointing vectors kx and ky

    :param mesh: [2,nx,ny] array of coordinates, first-axis corresponds to x and y axis,

    :param kx: horizontal pointing vector [-pi, pi]
    :param ky: vertical pointing vector [-pi, pi]
    """
    xx, yy = np.meshgrid(x, y)
    kx = 2 * np.pi * np.sin(theta_x) / e2wav(photon_energy)
    ky = 2 * np.pi * np.sin(theta_y) / e2wav(photon_energy)

    tilt = np.exp(1j * (kx * xx + ky * yy))
    tilt /= np.sqrt(np.sum(abs(tilt**2)))
    return tilt


def geometric_focus(fwhm, divergence):
    """
    calculate the focal length of a converging beam of size fwhm

    approximation that we are viewing the beam far downstream of the beam waist
    """
    return fwhm / np.tan(divergence)
