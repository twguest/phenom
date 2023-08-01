# -*- coding: utf-8 -*-

import numpy as np

from phenom.utils import e2wav


def gaussian_1d(x, a, mu, sigma):
    """
    generate a gaussian envelope for modelling the integrated pulse_data

    :param x: dependent-axis, one-dimensional numpy array
    :param a: amplitude of the gaussian (float)
    :param mu: mu
    :param width: widht of gaussian beam

    :returns: one-dimensioanl caussian profile
    """
    assert type(x) is np.ndarray

    return a * np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sigma, 2.0)))


def complex_gaussian_beam(x, y, photon_energy, pulse_energy, sigma, div, theta_x=0, theta_y=0, x0=0, y0=0):
    """
    Generates an analytical model for a complex Gaussian source, often used in optical simulations.

    :param x: 1d horizontal coordinate as a numpy array.
    :param y: 1d vertical coordinates as a numpy array.
    :param photon_energy: Photon energy in eV.
    :param pulse_energy: Pulse energy in Joules.
    :param sigma: Full width at half maximum (FWHM) of the Gaussian beam at a distance from the source.
    :param div: Beam divergence in radians.
    :param theta_x: Horizontal tilt wavevector (optional, default=0).
    :param theta_y: Vertical tilt wavevector (optional, default=0).
    :param x0: Horizontal shift (jitter) of the beam (optional, default=0).
    :param y0: Vertical shift (jitter) of the beam (optional, default=0).

    :returns: A complex Gaussian beam as a 2D numpy array.
    """
    xx, yy = np.meshgrid(x, y)

    R = (sigma) / div

    a = np.sqrt(pulse_energy)
    b = -((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma**2)
    c = ((xx - x0) ** 2 + (yy - y0) ** 2) / (4 * R)

    k = 2 * np.pi / e2wav(photon_energy)

    ### note, removed c and d (tilt)
    cgb = a * np.exp(b - 1j * k * c)
    cgb /= np.sqrt(np.sum(abs(cgb) ** 2))
    return cgb


if __name__ == "__main__":
    x = np.linspace(-5, 5, 100)
    g = gaussian_1d(x, a=1, mu=0, sigma=1)

    from matplotlib import pyplot as plt

    plt.plot(x, g)
    print(np.sum(g))
