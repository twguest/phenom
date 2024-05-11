# -*- coding: utf-8 -*-

import numpy as np

from matplotlib import pyplot as plt

from numpy import fft

from phenom.gaussian import gaussian_1d


def linear_SASE_spectrum(t, pulse_duration, photon_energy, bandwidth=1e-04, t0=0, norm=False, plot=False, output = 'freq'):
    """
    Generate a single SASE pulse profile in the time and frequency domain.

    :param t: Time axis [s](np.ndarray)
    :param pulse_duration: FWHM value of the SASE pulse duration [s](float)
    :param photon_energy: Central energy of the pulse [eV](float)
    :param bandwidth: Energy bandwidth / relative of the pulse [unitless](float)
    :param t0: Timing jitter (generally should be less than 2*sampling*bandwidth) [s](float)
    :param norm: Normalize intensity spectrum to 1 [default=False](bool)
    :param plot: Plot the outputs if True [default=False](bool)
    :param output: 'time' or 'frequency' spectrum 
    
    :returns: Complex spectrum in time or frequency domain [default time/s](np.ndarray)
    """
    hr_eV_s = 4.135667696e-15  # Planck constant in eV*s

    bw = photon_energy * bandwidth  # spectral width (FWHM)
    sigma_pulse = pulse_duration / (2 * np.sqrt(2 * np.log(2)))  # convert FWHM to standard deviation

    E = np.linspace(photon_energy - 3*bw, photon_energy + 3*bw, t.shape[-1])

    ### Gaussian Envelopes
    temporal_envelope = gaussian_1d(t, 1, t0, sigma_pulse)
    spectral_envelope = gaussian_1d(E, 1, photon_energy, bw)

    ### Random Phase and Fourier Transform
    random_phases = np.random.uniform(-np.pi, np.pi, E.size)
    complex_spectrum = spectral_envelope * np.exp(-1j * random_phases)
    ifft_spectrum = fft.ifftshift(fft.ifft(complex_spectrum))

    ### Applying Temporal Envelope
    time_domain_field = np.fft.ifftshift(ifft_spectrum * temporal_envelope)
    
    
    ### Fourier Transform back to Frequency Domain
    freq_domain_field = fft.fft(time_domain_field)  

    if plot:
        
        intensity_time_domain = abs(time_domain_field) ** 2
        intensity_freq_domain = abs(freq_domain_field) ** 2
    
        plt.figure(figsize=(10, 10))
        plt.subplot(411)
        plt.title("Temporal Envelope")
        plt.plot(t, temporal_envelope)
        plt.subplot(412)
        plt.title("Initial Intensity Spectrum")
        plt.plot(E, abs(complex_spectrum)**2)
        plt.subplot(413)
        plt.title("Time Domain Intensity")
        plt.plot(t, intensity_time_domain)
        plt.subplot(414)
        plt.title("Frequency Domain Intensity After Temporal Envelope")
        plt.plot(E, intensity_freq_domain)
        plt.tight_layout()
        plt.show()

    if output == 'time':
        complex_spectrum = time_domain_field
    elif output == 'freq':
        complex_spectrum = freq_domain_field
        
        
    return complex_spectrum


if __name__ == "__main__":
    t = np.linspace(-25e-15, 25e-15, 1500)
    spectrum = linear_SASE_spectrum(t, pulse_duration=5e-15, photon_energy=9500, bandwidth=1e-12, plot=True, output = 'time')
