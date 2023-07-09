# -*- coding: utf-8 -*-

import numpy as np
from numpy import fft
from scipy.constants import e, h

from phenom.gaussian import gaussian_1d

h_eV_s = h / e
hr_eV_s = h_eV_s / (2 * np.pi)

np.seterr(invalid="ignore")


def linear_SASE_spectrum(t, pulse_duration, photon_energy, bandwidth=1e-04, t0=0):
    """
    generate a single SASE pulse profiles

    assumes that the spectral and temporal profiles are gaussian,
    assumes that there is no jitter from the central value of the Gaussian (ie, the pulse profile
    is persistent).

    in future, we will extend this extned the functionality to account for non-Gaussian profiles,
    and pulse-length jitter.

    it will also be useful to trim the function to the relevant regions (and thus reduce the number of points)

    :param pulse_duration: expectation fwhm value of the SASE pulse time [in s]
    :param photon_energy: central energy of pulse [in eV]
    :param bandwidth: energy bandwidth / relative of the pulse
    :param t0: timing jitter (float) (should geberally be less than 2*sampling*bandwidth)
    :param t: time axis

    :returns: spectrum
    """
    nt = t.shape[0]

    pulse_duration = np.sqrt(2) * pulse_duration / (2 * np.sqrt(2 * np.log(2)))
    bw = photon_energy * bandwidth * np.sqrt(2)
    E = np.linspace(photon_energy - bw, photon_energy + bw, nt)  ### bandwidthfine frequency/energy domain

    estep = (E[1] - E[0]) / hr_eV_s  ### step size of freq domain

    et = np.linspace(-np.pi / estep, np.pi / estep, nt)  ### bandwidthfine time-domain

    temporal_envelope = (1 / np.sqrt(2 * np.pi)) * gaussian_1d(et, 1, t0, pulse_duration)

    spectral_envelope = gaussian_1d(E, 1, photon_energy, bw)

    random_phases = np.random.uniform(-np.pi, np.pi, temporal_envelope.shape)

    spectrum = fft.fftshift(fft.fft(spectral_envelope * np.exp(-1j * random_phases))) * temporal_envelope
    spectrum /= np.sqrt(np.sum(abs(spectrum)) ** 2)  ### normalise area under intensity curve to 1

    return spectrum


if __name__ == "__main__":
    spectrum = linear_SASE_spectrum(pulse_duration=25e-15, photon_energy=9500, bandwidth=1e-3, nt=750)
