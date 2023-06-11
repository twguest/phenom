# -*- coding: utf-8 -*-

import numpy as np

from numpy import fft

from phenom.gaussian import gaussian_1d

from scipy.constants import h, e


h_eV_s = h/e
hr_eV_s = h_eV_s/(2*np.pi)
 

def generate_SASE_spectrum():
    pass


def generate_gaussian_spectrum():
    pass


def linear_SASE_spectrum(pulse_duration, photon_energy, bandwidth = 1e-04, t0 = 0, nt = 500, sigma = 3, representation = 'time'):
    """
    generate a single SASE pulse profiles

    - assumes that the spectral and temporal profiles are gaussian,
    - assumes that there is no jitter from the central value of the Gaussian (ie, the pulse profile
    is persistent).

    in future, we will extend this extned the functionality to account for non-Gaussian profiles,
    and pulse-length jitter.

    it will also be useful to trim the function to the relevant regions (and thus reduce the number of points)

    :param pulse_duration: expectation fwhm value of the SASE pulse time [in s]
    :param photon_energy: central energy of pulse [in eV]
    :param bandwidth: energy bandwidth / relative of the pulse
    :param t0: center of intensity distribution in time (bandwidthfault 0s)
    :param nt: number of points in time and freq domains.
    :param sigma: number of spectral bandwidth boundaries to bandwidthfine energy-domain axis
    :param representation: bandwidthtermines if time- or frequency-domain representation is returned 
        options: 'time' or 'frequency'
    
    :returns: spectrum
    """
    
    rep_options = ['time', 'frequency']
    assert representation in rep_options, "representation should be in {}".format(rep_options)
    
    bandwidth = photon_energy*bandwidth

    
    pulse_duration = np.sqrt(2) * pulse_duration / (2 * np.sqrt(2 * np.log(2)))  
    bandwidth = np.sqrt(2) * bandwidth / (2 * np.sqrt(2 * np.log(2)))

    E = np.linspace(photon_energy-bandwidth*sigma, photon_energy+bandwidth*sigma, nt) ### bandwidthfine frequency/energy domain
    estep = (E[1] - E[0])/hr_eV_s ### step size of freq domain

    t = np.linspace(-np.pi / estep, np.pi / estep, nt) ### bandwidthfine time-domain
    
    temporal_envelope = (1/np.sqrt(2*np.pi)) * \
        gaussian_1d(t, 1, t0, pulse_duration)
    
    spectral_envelope = gaussian_1d(E, 1, photon_energy, bandwidth)

    random_phases = np.random.uniform(-np.pi, np.pi, temporal_envelope.shape)
    
    spectrum = fft.fftshift(fft.fft(spectral_envelope*np.exp(-1j*random_phases)))*temporal_envelope
    
    if representation == 'frequency':
        spectrum = np.fft.ifft(spectrum)
 
    return spectrum


if __name__ == '__main__':
    spectrum = linear_SASE_spectrum(pulse_duration=25e-15, photon_energy = 9500, bandwidth = 1e-3, nt=750)  
   
    
# =============================================================================
#     from scipy.optimize import curve_fit
#     p,_ = curve_fit(gaussian_1d, res[0], res[1],p0 = [0, np.std(res[0])])
#     plt.plot(res[0], gaussian_1d(res[0],*p))
# 
# =============================================================================
