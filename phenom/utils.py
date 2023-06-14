#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 20:03:10 2023

@author: twguest
"""

import numpy as np
from scipy.constants import c, e, h


def e2wav(photon_energy):
    return (h * c) / (e * photon_energy)


def e2k(photon_energy):
    return (np.pi * 2) / e2wav(photon_energy)
