#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:01:09 2023

@author: twguest
"""


from types import FunctionType

import h5py
import numpy as np

from phenom.gaussian import complex_gaussian_beam
from phenom.spectrum import linear_SASE_spectrum
from phenom.wavefront_tools import wavefront_tilt


def sase_pulse(
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
    Generate a Self-Amplified Spontaneous Emission (SASE) pulse model.

    This function models a SASE pulse by combining the temporal profile (from the linear_SASE_spectrum function), the spatial profile (from the complex_gaussian_beam function), and a wavefront tilt.

    :param x: numpy.ndarray
        Horizontal coordinates, representing the spatial dimension in the horizontal direction.

    :param y: numpy.ndarray
        Vertical coordinates, representing the spatial dimension in the vertical direction.

    :param t: numpy.ndarray
        Time coordinates, representing the temporal dimension of the pulse.

    :param photon_energy: float
        Photon energy of the pulse, in electronvolts (eV).

    :param pulse_energy: float
        Total energy of the pulse, in joules (J).

    :param pulse_duration: float
        Duration of the pulse, in seconds.

    :param bandwidth: float
        Bandwidth of the pulse, given as a fraction of the photon energy.

    :param sigma: float
        Gaussian beam width, representing the standard deviation of the spatial Gaussian profile.

    :param div: float
        Beam divergence, in radians.

    :param x0: float
        Horizontal position jitter, representing a shift in the horizontal direction.

    :param y0: float
        Vertical position jitter, representing a shift in the vertical direction.

    :param t0: float
        Temporal jitter, representing a shift in the temporal profile.

    :param theta_x: float
        Horizontal pointing angle, defining the tilt of the wavefront in the horizontal direction.

    :param theta_y: float
        Vertical pointing angle, defining the tilt of the wavefront in the vertical direction.

    :return: numpy.ndarray
        Electric field of the modeled SASE pulse, represented as a complex array.
    """

    tfield = linear_SASE_spectrum(
        t,
        pulse_duration=pulse_duration,
        photon_energy=photon_energy,
        bandwidth=bandwidth,
        t0=t0,
    )

    sfield = complex_gaussian_beam(
        x=x,
        y=y,
        photon_energy=photon_energy,
        pulse_energy=photon_energy,
        sigma=sigma,
        div=div,
        x0=x0,
        y0=x0,
        theta_x=theta_x,
        theta_y=theta_y,
    )

    tilt = wavefront_tilt(x=x, y=y, photon_energy=photon_energy, theta_x=theta_x, theta_y=theta_y)

    efield = sfield[:, :, np.newaxis] * tilt[:, :, np.newaxis] * tfield
    efield /= np.sqrt(np.sum(abs(efield) ** 2) * np.ptp(x) * np.ptp(y) * np.ptp(t))
    efield *= pulse_energy

    return efield


def check_arg_types(args):
    """
    Check the argument types in a given dictionary.
    
    This function verifies that each value in the dictionary corresponds to one of the allowed types. If the value is a list, it is converted to a NumPy array.
    
    :param args: dict
        A dictionary containing the arguments to be checked. The keys represent argument names, and the values are the corresponding objects to be verified.
    
    :raises AssertionError:
        If any value in the dictionary is not of one of the allowed types: float, list, numpy.ndarray, or FunctionType.
    
    :side effect:
        If any value in the dictionary is of type list, it is converted to a numpy.ndarray, and the dictionary is updated in place.
    
    :Note:
        This function is intended for internal use, to ensure that the inputs to a given function adhere to the expected types.
        It is useful in the development and debugging process, particularly when dealing with complex functions that accept multiple argument types.
    
    :Example:
        args = {
            'x': [1, 2, 3],
            'y': 5.0,
            'func': lambda x: x**2
        }
        check_arg_types(args)
        # args now contains: {'x': numpy.ndarray([1, 2, 3]), 'y': 5.0, 'func': <function>}
    """
    parse_types = [float, list, np.ndarray, FunctionType]

    keys = list(args.keys())

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
    Sorts the arguments based on their types and asserts array length consistency.

    This utility function organizes the arguments of the Source function by their types into different categories and ensures that all lists and 1D NumPy arrays have equal lengths. It's useful for categorizing and validating the input parameters.

    :param args: dict
        A dictionary containing the arguments of the Source function. The keys represent argument names, and the values are the corresponding objects to be sorted.

    :param __type__: type
        The expected types for the arguments.

    :returns N: int
        The number of pulses to iterate over. If there are no arrays in the arguments, this value will be 1. Otherwise, it will be the common length of all 1D NumPy arrays.

    :returns __set__: dict
        A dictionary with keys representing the types ("float", "np.ndarray", "FunctionType") and values as lists containing the corresponding keys from the input 'args' that match the type.

    :raises AssertionError:
        If the lengths of all lists and 1D NumPy arrays are not equal.

    :Example:
        args = {
            'x': np.array([1, 2, 3]),
            'y': 5.0,
            'func': lambda x: x**2
        }
        N, sorted_args = sort_arg_types(args, __type__)
        # N = 3
        # sorted_args = {'float': ['y'], 'np.ndarray': ['x'], 'FunctionType': ['func']}
    """
    __keys__ = list(args.keys())

    __arrays__ = []
    __lams__ = []

    N = 1  ### number of pulses

    for key in __keys__:
        # print(key,type(args[key])) # dev @author: twguest - 13/06/23

        if type(args[key]) == np.ndarray:
            __arrays__.append(key)

        elif type(args[key]) == FunctionType:
            __lams__.append(key)

    __len__ = [args[key].shape[0] for key in __arrays__]

    if len(__arrays__) > 0:
        assert np.all(__len__ == np.mean(__len__)), "all length of all lists and 1d numpy arrays must be equal"

        N = int(np.mean(__len__))

    __floats__ = list(set(__keys__) - set(__arrays__) - set(__lams__))

    __set__ = {}

    __set__["float"] = __floats__
    __set__["np.ndarray"] = __arrays__
    __set__["FunctionType"] = __lams__

    return N, __set__


def get_queue(args):
    """
    source models method of parsing input values to individual pulse values
    """
    try:
        del args["x"]
    except KeyError:
        pass

    try:
        del args["y"]
    except KeyError:
        pass

    try:
        del args["t"]
    except KeyError:
        pass

    __type__ = check_arg_types(args)
    N, __set__ = sort_arg_types(args, __type__)

    __queue__ = {}

    for key in __set__["float"]:
        __queue__[key] = np.ones([N]) * args[key]

    for key in __set__["FunctionType"]:
        assert type(args[key]()) == float, "all lambdas should generate a float"

        __queue__[key] = np.array([args[key]() for n in range(N)])

    for key in __set__["np.ndarray"]:
        __queue__[key] = args[key]

    return __queue__


def init_metadata():
    """
    wrapper to create metadata during Source init
    """
    metadata = {}
    metadata["stored"] = []


class Source:

    def __init__(self, **kwargs):
        """
        A generic class representing a source for generating electromagnetic wavefronts.

        Attributes:
        args : dict
            A dictionary containing the arguments that define the source properties.
        metadata : dict
            A dictionary containing metadata associated with the source.
        __queue__ : dict
            A dictionary representing the processing queue.
        N : int
            The number of processes in the queue.
        processes : dict
            A dictionary of processes representing individual pulses.

        Methods:
        get_process_parameters() -> None
            Calculates individual pulse parameters from the queue and stores them as processes.
        refresh_processes() -> None
            Regenerates the process list from stored arguments.
        execute(method, sdir) -> None
            A method to generate pulse wavefront files using a given method and save them to a specified directory.
        store_hdf5(sdir) -> None
            A method to store data in HDF5 format (Note: this method is incomplete in the provided code snippet).
        """
        self.args = locals()

        self.metadata = init_metadata()

        self.args.update(kwargs)

        del self.args["self"]
        del self.args["kwargs"]

        if "x" in list(self.args.keys()):
            self.x = self.args["x"]
        if "y" in list(self.args.keys()):
            self.y = self.args["y"]
        if "t" in list(self.args.keys()):
            self.t = self.args["t"]

        __queue__ = get_queue(self.args)

        for key in list(__queue__.keys()):
            assert type(__queue__[key]) == np.ndarray

        self.__queue__ = __queue__

        self.get_process_parameters()

    def get_process_parameters(self):
        """
        calculates set of individual pulse parameters from queue
        used to generate of processes from queue


        TO-DO: this function is where pulse name in the source .h5 is determined
        """

        if len(list(self.__queue__.keys())) == 0:
            self.N = 0
        elif len(list(self.__queue__.keys())) > 0:
            self.N = self.__queue__[list(self.__queue__.keys())[0]].shape[0]  ## get number of processes

        self.processes = {}

        for n in range(self.N):
            params = {}

            for key in list(self.__queue__.keys()):
                params[key] = self.__queue__[key][n]

            self.processes["pulse{:03}".format(n)] = params

        del self.__queue__

    def refresh_processes(self):
        """
        regeneerate processes list from stored args.

        used to regenerate pulses given input rules

        """
        if hasattr(self, "processes"):
            del self.processes

        __queue__ = get_queue(self.args)

        for key in list(__queue__.keys()):
            assert type(__queue__[key]) == np.ndarray

        self.__queue__ = __queue__
        self.get_process_parameters()

    def execute(self, method, sdir):
        """
        a non-parallelised method of generating pulse wavefront files using a method
        files are saved in perscribed sdirs

        :param method: wavefront generating method that takes in process params as its input
        :param sdir: single save directory

        ### execture needs a check that process keys match requirements of method

        """

        if type(sdir) == list:
            assert len(sdir) == self.N, "list of strings should be length N (number of pulses to be generated)"
        else:
            assert type(sdir) == str, "sdir should be type str or list"

        with h5py.File(sdir, "w") as hf:
            for itr, key in enumerate(list(self.processes.keys())):
                if type(sdir) == list:
                    file = sdir[itr]

                if type(sdir) == str:
                    if ".h5" in sdir:
                        file = sdir.split(".h5")[0] + "{:03}".format(itr) + ".h5"

                    if ".hdf5" in sdir:
                        file = sdir.split(".h5")[0] + "{:03}".format(itr) + ".hdf5"

                self.processes[key]["file"] = file

                efield = method(self.processes[key])

                group = hf.create_group(key)
                group.create_dataset("data", shape=efield.shape, data=efield)

                params = hf.create_group(key + "/params")

                for k in list(self.processes[key].keys()):
                    params.create_dataset(k, data=self.processes[key][k])

                mesh = hf.create_group(key + "/mesh")
                mesh.create_dataset("x", data=self.x)
                mesh.create_dataset("y", data=self.y)
                mesh.create_dataset("t", data=self.t)

            hf.close()
        ###

    def store_hdf5(self, sdir):
        self.metadata["stored"].append((sdir))

        hf = h5py.File(sdir, "w")


class SASE_Source(Source):
    def __init__(
        self,
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
        theta_y,
    ):
        """
        A derived class from Source, specifically representing a Self-Amplified Spontaneous Emission (SASE) source.
    
        Attributes:
        Inherits all attributes from the Source class.
    
        Methods:
        generate_sase_field(params) -> np.ndarray
            Generates a SASE field using the sase_pulse function with the specified parameters.
        generate_pulses(sdir) -> None
            Executes the generation of SASE fields for each process and saves the results in the specified directory.
    
        Parameters:
        x : np.ndarray
            Horizontal coordinates.
        y : np.ndarray
            Vertical coordinates.
        t : np.ndarray
            Time coordinates.
        photon_energy : float
            Photon energy in eV.
        pulse_energy : float
            Pulse energy in J.
        pulse_duration : float
            Pulse duration in seconds.
        bandwidth : float
            Pulse bandwidth as a fraction of photon energy.
        sigma : float
            Gaussian beam width.
        div : float
            Beam divergence.
        x0 : float
            Horizontal position jitter.
        y0 : float
            Vertical position jitter.
        t0 : float
            Temporal jitter.
        theta_x : float
            Horizontal pointing angle.
        theta_y : float
            Vertical pointing angle.
        """
        super().__init__(
            x=x,
            y=y,
            t=t,
            photon_energy=photon_energy,
            pulse_energy=pulse_energy,
            pulse_duration=pulse_duration,
            bandwidth=bandwidth,
            sigma=sigma,
            div=div,
            x0=x0,
            y0=y0,
            t0=t0,
            theta_x=theta_x,
            theta_y=theta_y,
        )

    def generate_sase_field(self, params):
        sase_field = sase_pulse(
            x=self.x,
            y=self.y,
            t=self.t,
            photon_energy=params["photon_energy"],
            pulse_energy=params["pulse_energy"],
            pulse_duration=params["pulse_duration"],
            bandwidth=params["bandwidth"],
            sigma=params["sigma"],
            div=params["div"],
            x0=params["x0"],
            y0=params["y0"],
            t0=params["t0"],
            theta_x=params["theta_x"],
            theta_y=params["theta_y"],
        )

        return sase_field

    def generate_pulses(self, sdir):
        self.execute(self.generate_sase_field, sdir)
