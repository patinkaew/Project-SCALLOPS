from MPA_component import *
from MPA_utils import *
from scipy.constants import Planck, speed_of_light

def pack_signal_pile(wavelength_pile, signal_matrix, input_area, dfdt):
    rslt = []
    for i, wl in enumerate(wavelength_pile):
        lst = []
        for complex_amp in signal_matrix[i]:
            amp, phase = rect_to_polar(complex_amp)
            energy = (amp ** 2) * dfdt
            lst.append(SignalPacket(wl, energy, input_area, dfdt, phase))
        rslt.append(lst)
    return np.transpose(np.array(rslt))

def unpack_signal_pile(signal_pile, dfdt):
    rslt = []
    for signal_t in signal_pile:
        lst = []
        for signal_packet in signal_t:
            amp = np.sqrt((signal_packet.E_in / dfdt))
            phase = signal_packet.phase
            lst.append(polar_to_rect(amp, phase))
        rslt.append(lst)
    return np.transpose(np.array(rslt))

def set_cross_section_parameter():
    # fit the cross section function
    # this is  for ti-sapphire
    xdata = np.array([600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025, 1050, 1075, 1100])
    ydata = np.array([0, 0.1, 0.3, 1 ,1.6, 2.8, 3.5, 3.8, 3.9, 3.8, 3.5, 3.0, 2.4, 1.9, 1.6, 1.2, 1, 0.8, 0.6, 0.5, 0.4])
    xdata = convert_unit(xdata, 'nm', 'm')
    ydata = convert_unit(ydata * (1e-19), 'cm^2', 'm^2')
    guess = [800e-9, 100e-9, 3.9e-23]
    return fit_emission_cross_section(xdata, ydata, guess)

def cross_section_func(wavelength, params, scaling = 1):
    return (scaling)*reciprocal_gaussian(wavelength, *params)

def calculate_crystal_J_sat(cross_section, freq):
    return (Planck * freq) / (cross_section_func(lambda2nu(freq)))

def plot_multi_pass_data(ax, passes, data_list, color_list, label_list, x_scaling = 1, y_scaling = 1, loc = ''):
    for path in range(passes):
        for idx, data in enumerate(data_list):
            x, y = data
            x = x_scaling * x
            y = y_scaling * y
            points_per_path = int(len(x) // passes)
            start = path*points_per_path
            end = (path+1)*points_per_path

            if not path == 0:
                ax.plot(x[start : end], y[start : end], c = color_list[idx])
            else: #first pass
                ax.plot(x[start : end], y[start : end], c = color_list[idx], label = label_list[idx])

            if path < passes - 1 and idx == len(data_list) - 1:
                ax.axvline(x = x[end] - 0.5*(x[1] - x[0]), ls = '--', c = 'black', alpha = 0.5)
    if loc == '':
        plt.legend()
    else:
        plt.legend(loc = loc)

def plot_multi_pass_last_value(ax, passes, data, color, format, x_scaling = 1, y_scaling = 1, x_offset = 0, y_offset = 0, fontsize = 10):
    x, y = data
    x = x_scaling * x
    y = y_scaling * y
    points_per_path = int(len(x) // passes)
    for path in range(passes):
        val = y[(path+1)*points_per_path - 1]
        ax.text(x[(path+1)*points_per_path - 1] + x_offset, val + y_offset, format.format(val), c = color, fontsize = fontsize)
