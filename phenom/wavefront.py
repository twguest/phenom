#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:01:52 2023

@author: twguest
"""

import numpy as np

from phenom.utils import e2k, e2wav

def wavefront_tilt(meshgrid, kx, ky):
    """ 
    define a two-dimensional complex tilted plane-wave, tilted along transverse pointing vectors kx and ky
    
    :param mesh: [2,nx,ny] array of coordinates, first-axis corresponds to x and y axis,
    see felpy.utils.np_utils.get_mesh for more
    :param kx: horizontal pointing vector [-pi, pi]
    :param ky: vertical pointing vector [-pi, pi]
    """
    
    return np.exp(1j*(kx*meshgrid[0]+ky*meshgrid[1]))


def spherical_phase(nx, ny,
                    xMin, xMax,
                    yMin, yMax,
                    fwhm, divergence,
                    photon_energy, x0=0, y0=0):

    k = e2k(photon_energy)
    f = geometric_focus(fwhm, divergence)

    # Initializing value of x-axis and y-axis
    # in the range -1 to 1
    x, y = np.meshgrid(np.linspace(xMin, xMax, nx),
                       np.linspace(yMin, yMax, ny))

    # Initializing sigma and muu

    # Calculating Gaussian array
    spherical_phase_term = np.exp(1j*k*((x-x0)**2+(y-y0)**2)/(4*f))

    return spherical_phase_term


def geometric_focus(fwhm, divergence):
    """
    calculate the focal length of a converging beam of size fwhm
    """
    return fwhm/np.tan(divergence)


class Wavefront:
    """
    wavefront is essentially a bare-bones version of wpg.wavefront.
    """
    
    def __init__(self, efield, nx, ny, nt, dx, dy, dt):
    
        self.mesh = {}
        self.mesh['nx'] = nx
        self.mesh['ny'] = ny
        self.mesh['nt'] = nt
        
        self.mesh['dx'] = dx
        self.mesh['dy'] = dy
        self.mesh['dt'] = dt
        
        self.efield = np.ones([nx,ny,nt])
        
 
    def get_axis(self, axis_name):
        
        if axis_name == 'x':
            axis = np.arange(-self.mesh['dx']*self.mesh['nx']//2,self.mesh['dx']*self.mesh['nx']//2,self.mesh['nx'])
        
        if axis_name == 'y':
            axis = np.arange(-self.mesh['dy']*self.mesh['ny']//2,self.mesh['dy']*self.mesh['ny']//2,self.mesh['ny'])

        
        if axis_name == 'z':
            axis = np.arange(-self.mesh['dt']*self.mesh['nt']//2,self.mesh['dt']*self.mesh['nt']//2,self.mesh['nt'])
        
        return axis
    
    def get_intensity(self):
        return abs(self.efield)**2
    
    
    
    

        
        