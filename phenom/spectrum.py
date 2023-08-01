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
    Generate a single SASE (Self-Amplified Spontaneous Emission) pulse profile.

    This function assumes that both the spectral and temporal profiles are Gaussian. The current implementation does not consider any jitter from the central value of the Gaussian, meaning the pulse profile is considered persistent.

    Future extensions to this function might include the ability to handle non-Gaussian profiles and pulse-length jitter, as well as the trimming of the function to relevant regions to reduce the number of points.

    :param t: numpy.ndarray
        Time axis of the pulse. Represents the time range over which the pulse is defined.

    :param pulse_duration: float
        Expected full width at half maximum (FWHM) value of the SASE pulse time, in seconds.

    :param photon_energy: float
        Central energy of the pulse, in electronvolts.

    :param bandwidth: float, optional
        Energy bandwidth relative to the pulse. Represents the spread in energy of the pulse.
        Default is 1e-04.

    :param t0: float, optional
        Timing jitter, which should generally be less than twice the product of the sampling rate and the bandwidth.
        Default is 0.

    :return: numpy.ndarray
        The complex spectrum of the SASE pulse, represented as a 1D array.

    :Description:
        The function begins by defining several parameters including the pulse duration in terms of standard deviation, the bandwidth in absolute terms, the energy range (E), and the corresponding time-domain (et).
        It then calculates the temporal and spectral envelopes using a Gaussian function (assumed to be `gaussian_1d` within the function's scope).
        Random phases are introduced to the spectral envelope, and the final spectrum is computed using a Fourier transform.
        The resulting spectrum is normalized so that the area under the intensity curve is 1.

    :Note:
        - The `gaussian_1d` function, used to define the temporal and spectral envelopes, should be available in the scope of this function.
        - Constants like `hr_eV_s` used for unit conversion should also be defined.
        - The provided time axis 't' should be a 1D array with an appropriate range to represent the pulse profile.
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


# =============================================================================
#     from scipy.optimize import curve_fit
#     p,_ = curve_fit(gaussian_1d, res[0], res[1],p0 = [0, np.std(res[0])])
#     plt.plot(res[0], gaussian_1d(res[0],*p))
#
# =============================================================================
