# -*- coding: utf-8 -*-

import numpy as np

from phenom.utils import e2wav

def gaussian_1d(x,a, mu, sigma):
    """
    generate a gaussian envelope for modelling the integrated pulse_data
    
    :param x: dependent-axis, one-dimensional numpy array
    :param a: amplitude of the gaussian (float)
    :param mu: mu
    :param width: widht of gaussian beam
    
    :returns: one-dimensioanl caussian profile
    """
    
    assert type(x) is np.ndarray 
    
    return a*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))


def complex_gaussian_beam(x,y,photon_energy,i,sig0,sigma,div,theta_x = 0, theta_y = 0,mu = 0, y0 = 0):
    """
    Analytical model for a complex gaussian source
    
    :param x: horizontal coordinates
    :param y: vertical coordinates
    :param i: on-axis beam intensity
    :param sig0: waist fwhm 
    :param sigma: fwhm @ a distance z from the source
    :param div: beam divergence
    :param theta_x: horizontal tilt wavevector
    :param theta_y: vertical tilt wavevector
    :param mu: horizontal shift 
    :param y0: vertical shift
    """

    R = (sigma*1.6)/div

    a = np.sqrt(i)*(sigma/sig0)
    b = -((x-mu)**2+(y-y0)**2)/(2*sigma**2)
    c = ((x-mu)**2+(y-y0)**2)/(4*R)
    
    k = (2*np.pi/e2wav(photon_energy))

    ### note, removed c and d (tilt)
    return a*np.exp(b-1j*k*c)




if __name__ == '__main__':
    x = np.linspace(-5,5,100)
    g = gaussian_1d(x, mu = 0, sigma = 1)
    
    from matplotlib import pyplot as plt
    plt.plot(x, g)
    print(np.sum(g))
    
    
    