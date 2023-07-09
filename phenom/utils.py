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
    Converts photon energy to photon wavelength.
    
    :param photon_energy: photon energy in eV
    
    :returns: photon wavelength in m
    """
    return (h * c) / (e * photon_energy)


def e2k(photon_energy):
    """
    Converts photon energy to freq. wavevector
    
    :param photon_energy: photon energy in eV
    
    :returns k: wavevector (1/m)
    """
    return (np.pi * 2) / e2wav(photon_energy)


def complex_to_wpg(arr): ### converter
    """
    converter function to transform complex wavefield into wpg style electric
    field array
    
    :param arr: complex wavefield array [x,y,t] (complex128 type)
    
    :returns new_arr: wpg style electric field array [nx,ny,nz,2] (float64)
    """
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr
