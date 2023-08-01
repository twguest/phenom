#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:03:10 2023

@author: twguest
"""

import numpy as np
from scipy.constants import c, e, h


def e2wav(photon_energy):
    """
    Convert photon energy to corresponding wavelength.
    
    :param photon_energy: float
        The photon energy in appropriate units (e.g., electronvolts).
    
    :return: float
        The corresponding wavelength in appropriate units (e.g., meters or angstroms).
    
    :Note:
        This function uses the Planck's constant (h), speed of light (c), and elementary charge (e)
        to perform the conversion. 
    """
    return (h * c) / (e * photon_energy)


def e2k(photon_energy):
    """
    Convert photon energy to corresponding wave number.
    
    :param photon_energy: float
        The photon energy in appropriate units (e.g., electronvolts).
    
    :return: float
        The corresponding wave number in appropriate units (e.g., rad/m).
    
    :Note:
        This function uses the `e2wav` function to first convert the photon energy to wavelength
        and then calculates the wave number as (2 * pi) / wavelength.
    """
    return (np.pi * 2) / e2wav(photon_energy)
