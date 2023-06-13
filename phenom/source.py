#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:01:09 2023

@author: twguest
"""

from copy import copy

import numpy as np
from types import FunctionType

import h5py_wrapper as h5w

from phenom.gaussian import complex_gaussian_beam
from phenom.wavefront_tools import wavefront_tilt, spherical_phase
from phenom.utils import e2wav

from phenom.spectrum import linear_SASE_spectrum




def sase_pulse(x,y,
               t,
               photon_energy,
               pulse_energy,
               pulse_duration,
               bandwidth,
               sigma,
               div,
               x0,
               y0,
               t0,
               theta_x,
               theta_y,
               ):
    
    """
    SASE pulse model
    
    :param x: horizontal coordinates (np.ndarray) 
    :param y: vertical coordinates (np.ndarray) 
    :param photon_energy: photon energy in eV (float)
    :param pulse_energy: pulse energy in J (float)
    :param pulse_duration: pulse duration in seconds (float)
    :param bandwidth: pulse bandwidth as a fraction of photon energy (float)
    :param sigma: gaussian beam width (float)
    :param div: beam divergence (float)
    :param x0: horizontal position jitter (float)
    :param y0: vertical position jitter (float)
    :param t0: temporal jitter (float)
    :param theta_x: horizontal pointing angle (float)
    :param theta_y: vertical pointing angle (float)
             
    """
    
    tfield = linear_SASE_spectrum(t,
                                  pulse_duration = pulse_duration,
                                  photon_energy = photon_energy,
                                  bandwidth = bandwidth,
                                  t0 = t0,
                                  )
    
    sfield = complex_gaussian_beam(x = x, y = y,
                                   photon_energy = photon_energy,
                                   pulse_energy = photon_energy,
                                   sigma = sigma,
                                   div = div,
                                   x0 = x0,
                                   y0 = x0,
                                   theta_x = theta_x,
                                   theta_y = theta_y)
    
    tilt = wavefront_tilt(x = x, y = y, 
                          photon_energy = photon_energy,
                          theta_x = theta_x,
                          theta_y = theta_y)
    
    efield = sfield[:,:,np.newaxis]*tilt[:,:,np.newaxis]*tfield
    efield/= np.sqrt(np.sum(abs(efield)**2) * np.ptp(x)*np.ptp(y)*np.ptp(t))
    efield*= pulse_energy

    return efield


def check_arg_types(args):
    
    parse_types = [float, list, np.ndarray, FunctionType]
    
    keys = (list(args.keys()))
    
    for key in keys:
        
        key_types = {}
        key_types[key] = type(args[key])
        
        # print(type(args[key])) # dev @author: twg - 13/06/23
        
        assert key_types[key] in parse_types, " input {} should be in {}".format(key, parse_types)

            
        ### list to array
        if type(args[key]) == list:
  
            args[key] = np.asarray(args[key])
            key_types[key] = type(args[key])
        
                        

                
def sort_arg_types(args, __type__):
    """ 
    parser utility: sorts function types + array length assertions
    
    :param args: arguments of the Source function
    :param __type__: argument types
    
    :returns N: number of pulses to iterate over
    :returns __set__: dict of argument keys sorted by value type
    """    

    __keys__ = list(args.keys())
    
    __arrays__ = []
    __lams__ = []
    
    N = 1 ### number of pulses
   
    for key in __keys__:
        
        
        # print(key,type(args[key])) # dev @author: twguest - 13/06/23
        
        if type(args[key]) == np.ndarray:
            
            __arrays__.append(key)
            
        elif type(args[key]) == FunctionType:
            __lams__.append(key)
                  
            
    __len__ = [args[key].shape[0] for key in __arrays__]
    
    if len(__arrays__) > 0:

        assert np.all( __len__ == np.mean(__len__)), "all length of all lists and 1d numpy arrays must be equal"
    
        N = int(np.mean(__len__))
        
    
    __floats__ = list(set(__keys__) - set(__arrays__) - set(__lams__))
    
    __set__ = {}
   
    __set__['float'] = __floats__
    __set__['np.ndarray'] = __arrays__
    __set__['FunctionType'] = __lams__
    
    return N, __set__
    
    
def get_queue(args):
    """ 
    source models method of parsing input values to individual pulse values
    """
    try:
        del args['x']
    except(KeyError):
        pass
    
    try:
        del args['y']
    except(KeyError):
        pass
    
    try:
        del args['t']
    except(KeyError):
        pass
        
                     
    __type__ = check_arg_types(args)
    N, __set__ = sort_arg_types(args, __type__)
    
    
    __queue__ = {}
    
    for key in __set__['float']:
        __queue__[key] = np.ones([N]) * args[key]
        
    for key in __set__['FunctionType']:

        
        assert type(args[key](0)) == float, "all lambdas should generate a float"
        
        __queue__[key] = np.array([args[key](n) for n in range(N)])

    return __queue__ 


class Source:
    """
    """
    

    def __init__(self, **kwargs):
 
        self.args = locals()
        
        self.args.update(kwargs)
        
        del self.args['self']
        del self.args['kwargs']
 
    
        if 'x' in list(self.args.keys()):
            self.x = self.args['x']
        if 'y' in list(self.args.keys()):
            self.y = self.args['y']
        if 't' in list(self.args.keys()):
            self.t = self.args['t']
                       
        __queue__ = get_queue(self.args)
        
            
        for key in list(__queue__.keys()):
            assert type(__queue__[key]) == np.ndarray
        
            
        self.__queue__ = __queue__


    def get_process_parameters(self):
        
        N = self.__queue__[list(self.__queue__.keys())[0]].shape[0] ## get number of processes
        
        self.processes = {}
        
        for n in range(N):
             
            params = {}
            
            for key in list(self.__queue__.keys()):
                params[key] = self.__queue__[key][n]
            
            if hasattr(self, 'x'):
                params['x'] = self.x
            if hasattr(self, 'y'):
                params['y'] = self.y
            if hasattr(self, 't'):
                params['t'] = self.t
                    
            self.processes[str(n)] = params

        del self.__queue__       
        
        
 
class SASE_Source(Source):

    def __init__(self,
                 x, 
                 y,
                 t,
                 photon_energy,
                 pulse_energy,
                 pulse_duration,
                 bandwidth,
                 sigma,
                 div,
                 x0,
                 y0,
                 t0,
                 theta_x,
                 theta_y):
        """
        initialisation function. 
        
        :param x: horizontal coordinates (np.ndarray) 
        :param y: vertical coordinates (np.ndarray) 
        :param photon_energy: photon energy in eV (float)
        :param pulse_energy: pulse energy in J (float)
        :param pulse_duration: pulse duration in seconds (float)
        :param bandwidth: pulse bandwidth as a fraction of photon energy (float)
        :param sigma: gaussian beam width (float)
        :param div: beam divergence (float)
        :param x0: horizontal position jitter (float)
        :param y0: vertical position jitter (float)
        :param t0: temporal jitter (float)
        :param theta_x: horizontal pointing angle (float)
        :param theta_y: vertical pointing angle (float)
        
        :param N: number of pulses - only valid if all other parameters are floats or lambdas
                      
        """
        super().__init__(x = x,
                         y = y,
                         t = t,
                         photon_energy = photon_energy,
                         pulse_energy = pulse_energy,
                         pulse_duration = pulse_duration,
                         bandwidth = bandwidth,
                         sigma = sigma,
                         div = div,
                         x0 = x0,
                         y0 = y0,
                         theta_x = theta_x,
                         theta_y = theta_y) 
    
        def generate_sase_field(params):
            
            sase_pulse(x = params['x'],
                             y = params['x'],
                             t = params['x'],
                             photon_energy = photon_energy,
                             pulse_energy = pulse_energy,
                             pulse_duration = pulse_duration,
                             bandwidth = bandwidth,
                             sigma = sigma,
                             div = div,
                             x0 = x0,
                             y0 = y0,
                             theta_x = theta_x,
                             theta_y = theta_y)
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
            parser
            
            s
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
 
                    
                    
    

