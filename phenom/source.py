#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:01:09 2023

@author: twguest
"""

import numpy as np


import h5py_wrapper as h5w

from phenom.gaussian import complex_gaussian_beam
from phenom.wavefront import wavefront_tilt, spherical_phase
from phenom.utils import e2wav

def generate_beam_envelope(self, pulse_properties):
        """
        a function to generate a gaussian photon field based on the source properties

        this can be extended at a later date to include perturbations, i.e the zernike polynomials.
        """

        xx, yy = np.meshgrid(np.linspace(pulse_properties['xMin'], pulse_properties['xMax'], pulse_properties['nx']), np.linspace(pulse_properties['yMin'], pulse_properties['yMax'], pulse_properties['ny']))
        
        return complex_gaussian_beam(x = xx, y = yy, photon_energy = pulse_properties['photon_energy'],
        i = pulse_properties[pulse_energy], sig0 = pulse_properties['sig0'],
        sigma = pulse_properties['sigma'], div = pulse_properties['div'],
        theta_x = pulse_properties['theta_x'], theta_y = pulse_properties['theta_y'],
        x0 = pulse_properties['x0'], y0 = pulse_properties['y0'])




def calc_pulse_energy(dx,dy,dt,E, pulse_intensity):
    """
    calculate energy of  in time domain

    :param wf: wpg.Wavefront structure
    :return: pulse energy value, J
    """
    J2eV = 6.24150934e18
        
    pulse_energy = pulse_intensity
    pulse_energy_J = pulse_energy * dx * dy * 1e6 * dt
    photons_per_pulse = pulse_energy_J * J2eV / E
    
    #print(
    #    "Number of photons per pulse: {:e}".format(photons_per_pulse)
    #)
    
    return pulse_energy_J, photons_per_pulse




def generate_pulse(temporal, spatial, pulse_energy):
    return spatial[:,:,np.newaxis]*temporal
    
# =============================================================================
#  
# class Source:
#     
#     def __init__(self, nx, ny, nt, dx, dy, dt):
#         """
#         """
# 
#         self.mesh = {}
#         self.mesh['nx'] = nx
#         self.mesh['ny'] = ny
#         self.mesh['nt'] = nt
#     
#         self.mesh['dx'] = dx
#         self.mesh['dy'] = dy
#         self.mesh['dz'] = dt
#         
#     
#     def generate():
#         """
#         this function generates a wavefront and assigns it certain properties
#         which in turn determine its statistics.
#         """
# =============================================================================
    
class Source:
    """
    a generalised source model
    """


    def __init__(self,
                mesh,
                photon_energy,
                pulse_energy,
                sig0,
                sigma,
                div,
                x0,
                y0,
                z0,
                pulse_duration,
                bandwidth,
                theta_x,
                theta_y,
                dR,
                spectrum = 'flat',
                spectrum_representation = 'time',
                mode = 0):

        """
        initialisation function. 
        
        :param photon_energy: photon energy in keV
        :param mesh: felpy style mesh file defining xMin, xMax, yMin, yMax, nx, ny
        :param x: horizontal coordinates
        :param y: vertical coordinates
        :param i: on-axis beam intensity
        :param sig0: waist fwhm 
        :param sigma: fwhm @ a distance z from the source
        :param div: beam divergence
        :param theta_x: horizontal tilt wavevector
        :param theta_y: vertical tilt wavevector
        :param x0: horizontal shift 
        :param y0: vertical shift
        :param pulse_duration: temporal duration of the pulse
        """
        
        self.metadata = {}
        self.metadata['pulses'] = []
        
        self.source_properties = {}
 
 
        self.mesh = mesh
    
        self.source_properties['photon_energy'] = photon_energy
        self.source_properties['pulse_energy'] = pulse_energy
        self.source_properties['sig0'] = sig0 
        self.source_properties['sigma'] = sigma
        self.source_properties['div'] = div
        self.source_properties['x0'] = x0
        self.source_properties['y0'] = y0
        self.source_properties['pulse_duration'] = pulse_duration
        self.source_properties['theta_x'] = theta_x
        self.source_properties['theta_y'] = theta_y
        self.source_properties['dR'] = dR
        self.source_properties['bandwidth'] = bandwidth
        self.source_properties['z0'] = z0
        
        self.spectrum = spectrum
        self.spectrum_representation  = spectrum_representation
        self.mode = mode
 
        
        for item in self.source_properties:
            assert type(self.source_properties[item]) is np.ndarray, "{} is a list - please input it as an np.ndarray instead"
    
        
        
    @property
    def _control(self):
        """
        a tool to return the longest array in the source properties (ie. the control variable)
        """
        d = {}
        for item in self.source_properties:
            if isinstance(self.source_properties[item], np.ndarray):
                d[item] = self.source_properties[item]
        try: 
            l = len(self.source_properties[max(d, key=lambda x:len(d[x]))])
        except(TypeError, ValueError):
            l = 1
        return l
    
    
    def build_properties(self):
        
        g = self._control
        
        for item in self.source_properties:
            
            if isinstance(self.source_properties[item], np.ndarray):

                if self.source_properties[item].shape[0] == g:
                    pass
                else:
                    self.source_properties[item] = np.repeat(self.source_properties[item], g)
                    
            elif type(self.source_properties[item]) is list:
                
                if len(self.source_properties[item]) == g:
                    pass
                else:
                    self.source_properties[item]*=np.ones(g).astype(type(self.source_properties[item][0]))
            
            else:
                self.source_properties[item] = np.repeat(self.source_properties[item], g)
                   
 
    def set_property(self, property_name, value):
        """
        wrapper function to store a source property
        
        :param property_name: name of the proeprty to be stored in dict 
        :type: str
        (
        :param value: value of the property
        :type: single value or list of values.
        """
        self.source_properties[property_name] = value


    def get_property(self, property_name):
        """
        return the value of a specific source property
    
        :param property_name: name of the proeprty to be read from dict 
        :type: str
        
        """
        assert property_name in self.source_properties, "Source property is not currently defined"
        return self.source_properties[property_name]
    
    
    def generate(self, array, pulse_properties, outdir = "./wfr_", save = True):

        
            
        wfr =  wavefront_from_array(array, nx=pulse_properties['nx'], ny=pulse_properties['ny'],
                                    nz=pulse_properties['nz'],
                                    dx=(pulse_properties['xMin']-pulse_properties['xMax'])/pulse_properties['nx'],
                                    dy=(pulse_properties['yMin']-pulse_properties['yMax'])/pulse_properties['ny'],
                                    dz=(pulse_properties['yMin']/pulse_properties['nz']),                                    
                                    photon_energy=pulse_properties['photon_energy'],
                                    pulse_duration=pulse_properties['pulse_duration'],
                                    source_properties = pulse_properties,
                                   elim = [pulse_properties['e_min'], pulse_properties['e_max']])
        
        
        wfr.scale_beam_energy(pulse_properties['pulse_energy'])
        
        self.metadata['pulses'].append(outdir)
        
        if save:
            wfr.store_hdf5(outdir)
        else:
            return wfr



    def store_hdf5(self, outdir, overwrite = False):
        """
        function to write source data( to a hdf5 file
        
        :param outdir: outdir of ".h5" file
        :type: str
        """
        
        if ".h5" or ".hdf5" in outdir:
            pass
        else: outdir+".h5"
        
        h5w.save(outdir, self.metadata, path = 'metadata/', write_mode = 'a', overwrite_dataset = overwrite)
        h5w.save(outdir, self.source_properties, path = 'source_properties/', write_mode = 'a', overwrite_dataset = overwrite)
        h5w.save(outdir, self.mesh, path = 'mesh/', write_mode = 'a', overwrite_dataset = overwrite)

    def load_hdf5(self, indir):
        """
        function to read data from a hdf5 file
        
        :param indir: directory to read ".h5" file from
        :type: str
        """
        self.metadata = h5w.load(indir,  path = 'metadata/')
        self.source_properties = h5w.load(indir,  path = 'source_properties/')
        self.mesh = h5w.load(indir,  path = 'mesh/')
       

    @property
    def pulses(self):
        """ 
        list of pulses generated by this source
        """
        return self.metadata['pulses']
    
    
    def generator(self, outdir = "", filename = None, N = 1, save = True):
            """
            this is the emission process, and generates N wavefronts according to the rules given by the source paramters file.
            
            currently, it will only support cases where the input dictionary has values in the form of a single value or list.
            
            :param
            """
            
            wavefronts = [] ### N*self._control length list of wavefronts
    
                
    
            for n in range(N):
    
                self.build_properties()
    
                for itr in range(self._control):
    
                    pulse_properties = {item:self.source_properties[item][itr] for item in self.source_properties}
                    t, tp, f, fp = self.get_temporal_profile(pulse_properties)
                    
 
                    tilt = wavefront_tilt(np.meshgrid(np.linspace(pulse_properties['xMin'], pulse_properties['xMax'], pulse_properties['nx']),
                                                    np.linspace(pulse_properties['yMin'], pulse_properties['yMax'], pulse_properties['ny'])),
                                        2*np.pi*np.sin(pulse_properties['theta_x'])/e2wav(pulse_properties['photon_energy']),
                                        2*np.pi*np.sin(pulse_properties['theta_y'])/e2wav(pulse_properties['photon_energy']))
                    
                    
                    spherical = spherical_phase(nx = pulse_properties['nx'],
                    ny= pulse_properties['ny'], 
                    xMin= pulse_properties['xMin'], 
                    xMax= pulse_properties['xMax'], 
                    yMin= pulse_properties['yMin'], 
                    yMax= pulse_properties['yMax'], 
                    fwhm= pulse_properties['sigma'], 
                    divergence= pulse_properties['div'], 
                    photon_energy= pulse_properties['photon_energy'])
                    
                    pulse_properties['e_min'] = f[0]
                    pulse_properties['e_max'] = f[-1]
                    
                    efield = self.generate_beam_envelope(pulse_properties)[:,:,np.newaxis]*tilt[:,:,np.newaxis]*fp
                    
                    if filename is None:
                        filename = "wfr_{}.h5".format(uuid.uuid4())
                                                  
                    if ".h5" or ".hdf5" in filename:
                        pass
                    else: filename+".h5"
                    
                    if save:
                        self.generate(efield, pulse_properties, outdir = outdir + filename, save = save)
                        filename = None
                    else:
                        wavefronts.append(self.generate(efield, pulse_properties, outdir = outdir + filename, save = save))
                        filename = None
                        
            if save is False:
                return wavefronts
        
        
       
        


    def generate_beam_envelope(self, pulse_properties):
        """
        a function to generate a gaussian photon field based on the source properties

        this can be extended at a later date to include perturbations, i.e the zernike polynomials.
        """

        xx, yy = np.meshgrid(np.linspace(pulse_properties['xMin'], pulse_properties['xMax'], pulse_properties['nx']), np.linspace(pulse_properties['yMin'], pulse_properties['yMax'], pulse_properties['ny']))
        
        return complex_gaussian_beam(x = xx, y = yy, photon_energy = pulse_properties['photon_energy'],
        i = pulse_properties[pulse_energy], sig0 = pulse_properties['sig0'],
        sigma = pulse_properties['sigma'], div = pulse_properties['div'],
        theta_x = pulse_properties['theta_x'], theta_y = pulse_properties['theta_y'],
        x0 = pulse_properties['x0'], y0 = pulse_properties['y0'])





    def generator(self, outdir = "", filename = None, N = 1, save = True):
        """
        this is the emission process, and generates N wavefronts according to the rules given by the source paramters file.
        
        currently, it will only support cases where the input dictionary has values in the form of a single value or list.
        
        :param
        """
        
        wavefronts = [] ### N*self._control length list of wavefronts

            

        for n in range(N):

            self.build_properties()

            for itr in range(self._control):

                pulse_properties = {item:self.source_properties[item][itr] for item in self.source_properties}
                t, tp, f, fp = self.get_temporal_profile(pulse_properties)
                
                if 'nz' in pulse_properties:
                    pass
                else:
                    pulse_properties['nz'] = len(tp)

                tilt = wavefront_tilt(np.meshgrid(np.linspace(pulse_properties['xMin'], pulse_properties['xMax'], pulse_properties['nx']),
                                                np.linspace(pulse_properties['yMin'], pulse_properties['yMax'], pulse_properties['ny'])),
                                    2*np.pi*np.sin(pulse_properties['theta_x'])/photon_energy2wav(pulse_properties['photon_energy']),
                                    2*np.pi*np.sin(pulse_properties['theta_y'])/photon_energy2wav(pulse_properties['photon_energy']))
                
                
                spherical = spherical_phase(nx = pulse_properties['nx'],
                ny= pulse_properties['ny'], 
                xMin= pulse_properties['xMin'], 
                xMax= pulse_properties['xMax'], 
                yMin= pulse_properties['yMin'], 
                yMax= pulse_properties['yMax'], 
                fwhm= pulse_properties['sigma'], 
                divergence= pulse_properties['div'], 
                photon_energy= pulse_properties['photon_energy'])
                
                pulse_properties['e_min'] = f[0]
                pulse_properties['e_max'] = f[-1]
                
                efield = self.generate_beam_envelope(pulse_properties)[:,:,np.newaxis]*tilt[:,:,np.newaxis]*fp
                
                if filename is None:
                    filename = "wfr_{}.h5".format(uuid.uuid4())
                                              
                if ".h5" or ".hdf5" in filename:
                    pass
                else: filename+".h5"
                
                if save:
                    self.generate(efield, pulse_properties, outdir = outdir + filename, save = save)
                    filename = None
                else:
                    wavefronts.append(self.generate(efield, pulse_properties, outdir = outdir + filename, save = save))
                    filename = None
                    
        if save is False:
            return wavefronts

            
 
        
    def get_temporal_profile(self, pulse_properties, sigma=4, refresh=False):
        ### note this gets and sets - should be changed
    
        if 'nz' in pulse_properties:
            n_samples = pulse_properties.get('nz')
        else:
            n_samples = 150
            n_samples *= sigma
        
        self.source_properties['nz'] = n_samples * np.ones(self._control).astype(type(n_samples))
        
        t, ET, f, EF = linear_SASE_spectrum(pulse_duration = pulse_properties['pulse_duration'],
                                    E0 = pulse_properties['photon_energy']*1e3, dE = pulse_properties['bandwidth'], sigma = 3, n_samples = n_samples)
        return t, ET, f, EF




class SA1_Source(Source):
    """
    Analytical SA1 Source Model
    """    

    
    def __init__(self, photon_energy, q, nx = 512, ny = 512, theta_x = 0, theta_y = 0, x0 = 0, y0 = 0, z0 = 0, mode = 0, bandwidth = 1e-03, **kwargs):

        """
        initialisation function. 

        :param photon_energy: photon energy in keV [float or list of floats of length N]
        :param q: beam charge in nC [float or list of floats of length N]
        :param nx: number of horizontal grid/mesh points (int)
        :param ny: number of vertical grid/mesh points (int)
        :param theta_x: horizontal pointing error [float or list of floats of length N]
        :param theta_y: vertical pointing error [float or list of floats of length N]
        :param x0: horizontal beam center [float or list of floats of length N]
        :param y0: vertical beam center [float or list of floats of length N]
        :param z0: position from source waist (float64)
        :param mode: defines how multi-dimensional inputs are treated
        :param bandwidth: energy bandwidth (fwhm of E**2) in eV
        
        if mode 0: arrays are taken to represent values of successive pulses
        if mode 1: arrays are taken to represent parameter space to scan over and a meshgrid is created.
            mode 1 is currently only supported for theta_x, theta_y

        :keyword i: pulse energy (J) [float or list of floats of length N]
        :keyword div: beam divergence (rad) [float or list of floats of length N]
        :keyword sig0: beam size at waist (m) [float or list of floats of length N]
        """
        
        if mode == 1:
            theta_x, theta_y, z = np.meshgrid(theta_x, theta_y, z0)
            theta_x = theta_x.reshape(functools.reduce(operator.mul, theta_x.shape, 1))
            theta_y = theta_y.reshape(functools.reduce(operator.mul, theta_y.shape, 1))
        

        self.source_properties = {}
        
        self.metadata = {}
        self.metadata['pulses'] = []

        self.source_properties['photon_energy'] = photon_energy
        self.source_properties['z0'] = z0

        if pulse_energy in kwargs:
            self.source_properties[pulse_energy] =  kwargs.get(pulse_energy)
            self.source_properties['pulse_energy'] = kwargs.get(pulse_energy)
        else:
            self.source_properties[pulse_energy] = analytical_pulse_energy(q, photon_energy) 
            self.source_properties['pulse_energy'] = analytical_pulse_energy(q, photon_energy) 

        if 'nz' in kwargs:
            self.source_properties['nz'] = kwargs.get('nz')
        
        if 'sig0' in kwargs:
            self.source_properties['sig0'] = kwargs.get('sig0')
        else:
            self.source_properties['sig0'] = analytical_pulse_width(photon_energy) 
        
        if 'div' in kwargs:
            self.source_properties['div'] = kwargs.get('div')
        else:
            self.source_properties['div'] = analytical_pulse_divergence(photon_energy, 'mean')
        
        self.source_properties['bandwidth'] = bandwidth

        self.source_properties['sigma'] =  self.source_properties['sig0'] + (z0*analytical_pulse_divergence(photon_energy, 'mean')*2)
        
        self.source_properties['fwhm (source)'] = self.source_properties['sig0']   
        self.source_properties['fwhm (z)'] = self.source_properties['sigma']   

        self.source_properties['sigma'] /= FWHM2E2
        self.source_properties['sig0'] /= FWHM2E2   
        
        if theta_x is list:
            self.source_properties['theta_x'] = [-t for t in theta_x]
        else:
            self.source_properties['theta_x'] = theta_x
        
        if theta_y is list:
            self.source_properties['theta_y'] = [-t for t in theta_y]
        else:
            self.source_properties['theta_y'] = theta_y   

        ### geometry approximation:
        self.source_properties['x0'] = x0 + np.tan(theta_x)*z0
        self.source_properties['y0'] = y0 + np.tan(theta_y)*z0
        
        

        if 'pulse_duration' in kwargs:
            self.source_properties['pulse_duration'] = kwargs.get('pulse_duration')
        else:
            self.source_properties['pulse_duration'] = analytical_pulse_duration(q)
 
        for item in self.source_properties:

            if type(self.source_properties[item]) == list:
                self.source_properties[item] = np.asarray(self.source_properties[item])
        
        ### approximations:
        # let -xMin = xMax = 1.5w, likewise for vertical coordinates
        r = 5
        self.source_properties['xMin'] = r*FWHM2E2*self.source_properties['sigma']+self.source_properties['x0'] 
        self.source_properties['yMin'] = r*FWHM2E2*self.source_properties['sigma']+self.source_properties['y0']

        self.source_properties['xMax'] = -r*FWHM2E2*self.source_properties['sigma']+self.source_properties['x0']
        self.source_properties['yMax'] = -r*FWHM2E2*self.source_properties['sigma']+self.source_properties['y0']
        
        self.source_properties['nx'] = nx
        self.source_properties['ny'] = ny
        
    @property
    def photon_energy(self):
        """
        return photon beam energy
        """
        return self.source_properties['photon_energy']
    
    
    @property
    def q(self):
        """
        return photon beam charge
        """
        return self.source_properties['q']
 
                    
                    
    

