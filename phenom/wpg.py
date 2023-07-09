#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:10:23 2023

All functions in this file require the install of WPG

@author: twguest
"""
import numpy as np
import h5py as h5

from wpg.srw import srwlpy 

from wpg.wavefront import Wavefront 

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


def wpg_converter(save_loc, key):
    """
    Convert saved phenom pulse to wpg type pulse
    
    :param save_loc: location of source .h5 master file (str)
    :param key: pulse key - e.g. 'pulse000' (str)
    
    :returns wfr: wpg compatible wavefront file
    """
    with h5.File(save_loc, mode = 'r') as hf:

        efield = hf[key]['data'][()]
        x = hf[key]['mesh']['x'][()]
        y = hf[key]['mesh']['y'][()]
        t = hf[key]['mesh']['t'][()]
        photon_energy = hf['pulse000']['params']['photon_energy'][()]
    
    
    nx, ny, nt = efield.shape

    wfr = Wavefront()


    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nt, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nt, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = photon_energy
    wfr.params.wDomain = 'time'
    wfr.params.Mesh.nSlices = nt
    wfr.params.Mesh.nx = nx
    wfr.params.Mesh.ny = ny

    wfr.params.Mesh.sliceMin = np.min(t)
    wfr.params.Mesh.sliceMax = np.max(t)


    wfr.params.Mesh.xMin = np.min(x)
    wfr.params.Mesh.xMax = np.max(x)
    wfr.params.Mesh.yMin = np.min(y)
    wfr.params.Mesh.yMax = np.max(y)

    wfr.params.Rx = 1
    wfr.params.Ry = 1

    wfr.data.arrEhor = complex_to_wpg(efield)
    srwlpy.SetRepresElecField(wfr._srwl_wf, 'frequency')
    return wfr