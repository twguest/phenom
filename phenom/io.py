#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:10:23 2023

@author: twguest
"""
import h5py as h5


def load_hdf5(file, proc):
    """
    function to read data from a hdf5 file to create a dictionary for each pulse in a source file

    :param file: source .h5 directory (str)
    :param proc: sub-group within .h5 (str)
    """
    wf = {}

    wf["mesh"] = {}
    wf["params"] = {}

    f = h5.File("./tmp/sase_field.h5", "r")

    assert proc in list(f.keys()), "sub-group '{}' does not exist".format(proc)

    wf["efield"] = f[proc]["data"]

    for key in list(f[proc]["params"].keys()):
        wf["params"][key] = f[proc]["params"][key]

    for key in list(f[proc]["mesh"].keys()):
        wf["mesh"][key] = f[proc]["mesh"][key]

    return wf
