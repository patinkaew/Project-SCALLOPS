import numpy as np
from scipy import optimize, fftpack, signal
from scipy.constants import nu2lambda, lambda2nu
import matplotlib.pyplot as plt
from pint import UnitRegistry

def polar_to_rect(spectral_amplitude, spectral_phase):
    return np.multiply(spectral_amplitude, np.exp(1j*spectral_phase))

def rect_to_polar(complex_amplitude):
    amp = np.absolute(complex_amplitude)
    phase = np.angle(complex_amplitude)
    return amp, phase

def gaussian(x, mean, sigma, amp = 1):
    return amp * (np.e ** ((-1/2) * ((x - mean) / sigma)**2))

def reciprocal_gaussian(x, mean, sigma, amp = 1):
    return amp * (np.e ** ((-1/2) * ((1/x - 1/mean) / (sigma/(mean * mean)))**2))

def fit_emission_cross_section(xdata, ydata, guess):
    params, cov = optimize.curve_fit(reciprocal_gaussian, xdata, ydata, p0 = guess)
    return params

def convert_unit(value, from_unit, to_unit):
    ureg = UnitRegistry()
    quantity = value * ureg(from_unit)
    return quantity.to(ureg(to_unit)).magnitude

def amplitude_to_energy(amp, df, dt):
    return np.absolute(amp) * np.absolute(amp) * df * dt #energy = |amplitude|^2 df dt

def energy_to_amplitude(energy, df, dt):
    return np.sqrt(energy/(df * dt))

def find_nearest(arr, val):
    index = (np.abs(np.array(arr) - val)).argmin()
    return arr[index]

def find_nearest_arg(arr, val):
    return (np.abs(np.array(arr) - val)).argmin()

def find_FWHM_arg(x, arr): # assume that the maximum is in the middle
    m_idx = np.argmax(arr)
    half_max = arr[m_idx]/2
    l_idx = find_nearest_arg(arr[:m_idx], half_max)
    r_idx = m_idx + find_nearest_arg(arr[m_idx:], half_max)
    return x[l_idx], x[r_idx]

def find_FWHM(x, arr): # assume that the maximum is in the middle
    m_idx = np.argmax(arr)
    half_max = arr[m_idx]/2
    l_idx = find_nearest_arg(arr[:m_idx], half_max)
    r_idx = m_idx + find_nearest_arg(arr[m_idx:], half_max)
    return np.abs(x[r_idx] - x[l_idx])
