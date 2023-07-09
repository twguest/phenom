#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:10:23 2023

@author: twguest
"""
import h5py as h5


def load_hdf5(file, pulse_id):
    """
    function to read data from a hdf5 file to create a dictionary for each pulse in a source file

    :param file: source .h5 directory (str)
    :param pulse_id: sub-group within .h5 (str)
    """
    wf = {}

    wf["mesh"] = {}
    wf["params"] = {}

    f = h5.File(file, "r")

    assert pulse_id in list(f.keys()), "sub-group '{}' does not exist".format(pulse_id)

    wf["efield"] = f[pulse_id]["data"]

    for key in list(f[pulse_id]["params"].keys()):
        wf["params"][key] = f[pulse_id]["params"][key]

    for key in list(f[pulse_id]["mesh"].keys()):
        wf["mesh"][key] = f[pulse_id]["mesh"][key]

    f.close()
    return wf


