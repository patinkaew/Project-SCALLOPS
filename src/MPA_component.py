from MPA_utils import *
import numpy as np
from scipy.constants import pi, speed_of_light, nu2lambda
from scipy import fftpack, interpolate, signal
from pandas import read_csv, DataFrame
import random

class Crystal():
    ''' Crystal class contains information about crystal used in amplification '''

    # fields of Crystal class
    __slots__ = ['attenuation', 'qd', 'length', 'n2', 'J_sat']

    # constructor
    def __init__(self, attenuation, qd, length, n2, J_sat):
        self.attenuation = attenuation # attenuation 'alpha' in Beer-Lambert law
        self.qd = qd # quantum defect
        self.length = length # length of the crystal
        self.n2 = n2 # nonlinear refractive index
        self.J_sat = J_sat # saturation fluence

    # create a copy of crystal with the same information
    def copy(self):
        return Crystal(self.attenuation, self.qd, self.length, self.n2, self.J_sat)

class Pump():
    # fields of Pump class
    __slots__ = ['E_pump', 'A_pump', 'J_pump']

    # constructor
    def __init__(self, E_pump, A_pump):
        self.E_pump = E_pump # pump energy
        self.A_pump = A_pump # pump area
        self.J_pump = E_pump / A_pump # pump fluence

class Store():
    # fields of Store class
    __slots__ = ['E_sto', 'A_sto', 'J_sto']

    # constructor
    def __init__(self, E_sto, A_sto):
        self.E_sto = E_sto # stored energy
        self.A_sto = A_sto # store area
        self.J_sto = E_sto / A_sto # stored fluence

class SignalPacket():
    __slots__ = ['wavelength', 'E_in', 'A_in', 'J_in', 'dfdt', 'phase', 'B_integral']
    def __init__(self, wavelength, E_in, A_in, dfdt, phase, B_integral = 0):
        self.wavelength = wavelength # wavelength of signal packet
        self.E_in = E_in # energy in the signal packet
        self.A_in = A_in # area of signal
        self.J_in = E_in / A_in # input fluence = input energy per area
        self.dfdt = dfdt # for converting back to intensity
        self.phase = phase # phase for dispersion calculation
        self.B_integral = B_integral # B integral

    # overloading = operation
    def __eq__(self, other):
        return self.__slots__ == other.__slots__

    # create a copy of signal packet with the same information
    def copy(self):
        return SignalPacket(self.wavelength, self.E_in, self.A_in, self.dfdt, self.phase, self.B_integral)

class SignalGenerator():
    def __init__(self, signal_res = 100):
        self.signal_res = int(2 * (signal_res //2))

    def set_signal_res(self, res):
        if res >= 1: #if the given resolution is invalid, this is unchanged
            self.signal_res = res
        return self.signal_res

    def calculate_pulse_energy(self, time, signal):
        return np.trapz(signal, dx = time[1] - time[0])

    def calculate_peak_power(self, pulse_duration, pulse_energy, signal_type, **kwargs): #from given pulse duration and pulse power
        time, signal = self.generate_signal(pulse_duration, 1, signal_type, **kwargs)
        basic = np.trapz(signal, dx = time[1] - time[0])
        return pulse_energy / basic

    def generate_signal(self, pulse_duration, peak_power, signal_type, pad_zero = True, xoffset = 0, yoffset = 0, **kwargs):
        time = np.linspace(-pulse_duration/2, pulse_duration/2, self.signal_res) - xoffset

        if signal_type == 'gauss':
            if 'FWHM' in kwargs:
                pulse = 2 **((-1/2) * (time / (kwargs['FWHM'] / 2))**2)
            else:
                raise Exception('missing FWHM for gaussian wave')
        elif signal_type == 'sech2':
            if 'FWHM' in kwargs:
                pulse = 1 / (np.cosh(time / (1.76 * kwargs['FWHM'] / 2)) ** 2)
            else:
                raise Exception('missing FWHM for sech square wave')
        elif signal_type == 'square' or signal_type == 'sq':
            pulse = np.ones(self.signal_res)
        elif signal_type == 'sawtooth' or signal_type == 'saw':
            pulse = 1/pulse_duration * (time + xoffset) + 0.5
        elif signal_type == 'triangle' or signal_type == 'tri':
            pulse = 1 - np.abs(2/pulse_duration * (time + xoffset))
        elif signal_type == 'reverse sawtooth' or  signal_type == 'revsaw':
            pulse = -1/pulse_duration * (time + xoffset) + 0.5

        elif signal_type == 'interpolate':
            if 'xpoints' in kwargs and 'ypoints' in kwargs and 'kind' in kwargs:
                xp = kwargs['xpoints']
                max = np.amax(xp)
                min = np.amin(xp)
                mid_range = (max + min) / 2
                range = max - min
                xpoints = (pulse_duration / range) * (xp - mid_range)
                f = interpolate.interp1d(xpoints, kwargs['ypoints'], kind = kwargs['kind'])
                pulse = f(time)
                pulse = pulse / np.max(pulse) #scaling
            else:
                raise Exception('missing interpolation information')

        elif signal_type == 'random':
            if 'random_range' in kwargs:
                min, max  = kwargs['random_range']
                if min < 0 or max > 1 :
                    raise ValueError('range must be between 0 and 1')
            else:
                min = 0
                max = 1
            pulse = np.array([min + random.random() * (max - min) for _ in range(self.signal_res)])
        else:
            raise ValueError('Incorrect signal type')

        #scaling with the peak power
        pulse = peak_power * (pulse + yoffset) / (1 + yoffset)

        #pad zero before and after the pulse
        if pad_zero:
            zero_pad = np.zeros(self.signal_res // 2)
            rslt = np.append(zero_pad, pulse)
            rslt = np.append(rslt, zero_pad)
            return np.linspace(-pulse_duration, pulse_duration, len(rslt)), rslt
        else:
            return time, pulse

class GratingPair():
    # fields of GratingPair class
    __slots__ = ['grating_L', 'grating_d']

    # constructor
    def __init__(self, grating_separation, grating_constant):
        self.grating_L = grating_separation # separation between two gratings
        self.grating_d = grating_constant # groove spacing of grating

    # add phase dispersion by grating pair dispersion equation
    def add_dispersion(self, incident_angle, input_time, input_complex_time, center_freq, mixing = 1):
        theta_i = np.deg2rad(incident_angle) # convert angle in degree to radian

        d_time = np.abs(input_time[1] - input_time[0]) # time spacing
        d_freq = 1/d_time # frequency spacing
        numpoints = len(input_complex_time) # number of point in the input signal
        input_complex_freq = fftpack.fft(input_complex_time) # fourier transform to get frequency domain
        input_freq = fftpack.fftfreq(numpoints) * d_freq # calculate frequency bins
        input_freq = np.fft.fftshift(input_freq) # rearrange frequency bins
        input_freq += center_freq # center frequency bins at center frequency
        input_amp_freq, input_phase_freq = rect_to_polar(input_complex_freq) # split into amplitude and phase

        center_wavelength = nu2lambda(center_freq) # wavelength corresponding to center frequency
        theta_r0 = np.arcsin(center_wavelength/self.grating_d - np.sin(theta_i)) # reflected angle of center wavelength
        k_0 = 2 * np.pi/center_wavelength * self.grating_L # phase velocity shift
        reciprocal_vg = ((1 + np.sin(theta_i) * np.sin(theta_r0))/(speed_of_light * np.cos(theta_r0)))*self.grating_L # group velocity shift

        # add phase to input signal using grating dispersion equation
        add_phase = [] # phase to be added
        output_amp_freq = input_amp_freq # output amplitude in frequency domain
        for index, freq in enumerate(input_freq):
            wavelength = nu2lambda(freq)
            dw = 2*np.pi*(freq - center_freq)
            sin_theta_r = wavelength / self.grating_d - np.sin(theta_i) # sine of reflected angle
            if np.abs(sin_theta_r) <= 1: # value of sine is valid
                add_ph = 2 * np.pi / (wavelength) * self.grating_L * np.cos(np.arcsin(sin_theta_r)) # grating-added dispersion
                add_ph -= (k_0 + reciprocal_vg*dw) # remove pulse shifting (only want changes of pulse's shape)
            else: # value is invalid
                add_ph = 0
                output_amp_freq[index] = 0 # this wavelength is not reflected
            add_phase.append(add_ph)
        add_phase = np.array(add_phase)

        output_phase_freq = input_phase_freq + mixing * add_phase # add phase to the input pulse

        output_complex_freq = polar_to_rect(output_amp_freq, output_phase_freq)
        output_complex_time = fftpack.ifft(output_complex_freq) # inverse fourier transform to get time domain

        return output_complex_time

    # stretch pulse one time (mixing = 1)
    def stretch_pulse(self, incident_angle, input_time, input_complex_time, center_freq):
        return self.add_dispersion(incident_angle, input_time, input_complex_time, center_freq, 1)

    # compress pulse one time (mixing = -1)
    def compress_pulse(self, incident_angle, input_time, input_complex_time, center_freq):
        return self.add_dispersion(incident_angle, input_time, input_complex_time, center_freq, -1)

class SignalManipulator():
    def __init__(self):
        pass

    def manipulate_signal_pile(self, signal_pile, energy_scaling, energy_offset, area_scaling, area_offset):
        new_signal_pile = []
        for signal_stack in signal_pile:
            new_signal_stack = []
            for signal_packet in signal_stack:
                new_signal_packet = signal_packet.copy()
                new_signal_packet.E_in =  energy_scaling * new_signal_packet.E_in + energy_offset
                new_signal_packet.A_in =  area_scaling * new_signal_packet.A_in + area_offset
                new_signal_stack.append(new_signal_packet)
            new_signal_pile.append(new_signal_stack)
        return new_signal_pile

    def divide_energy_signal_pile(self, signal_pile, denominator):
        return self.manipulate_signal_pile(signal_pile, 1/denominator, 0, 1, 0)

    def multiply_energy_signal_pile(self, signal_pile, multiplier):
        return self.manipulate_signal_pile(signal_pile, multiplier , 0, 1, 0)

class UnitAmplifier():
    # fielrs of unit amplifier
    __slots__ = ['crystal', 'store', 'input_stack', 'output_stack', 'E_rem']

    # contructor
    def __init__(self, crystal:Crystal, store:Store, signal_stack):
        self.crystal = crystal
        self.store = store

        if type(signal_stack) is SignalPacket: # only one packet
            signal_stack = [signal_stack]

        if not type(signal_stack[0]) is SignalPacket: #check error
            raise TypeError('Expect signel stack, a list of SignalPacket')
        self.input_stack = signal_stack

    def calculate_output(self):
        self.E_rem = self.store.E_sto
        output_stack = []
        for packet in self.input_stack:
            # calculate saturation fluence to account for gain narrowing
            J_sat = 0
            if type(self.crystal.J_sat) is dict or type(self.crystal.J_sat) is np.ndarray: # for multi wavelength
                J_sat = self.crystal.J_sat[packet.wavelength]
            else: # for single wavelength
                J_sat = self.crystal.J_sat

            # calculate output fluence
            if J_sat == np.infty or J_sat == 0:
                # J_sat is infinite if cross section is zero (far from center wavelength)
                # J_sat is zero if frequency of signal is zero (uncommon)
                J_out = packet.J_in # gain is one in both cases
            else:
                #calculate gain coefficient
                G_0 = np.exp(self.store.J_sto / J_sat)
                #calculate output fluence using Frantz - Nodvick equation
                J_out = J_sat * np.log(1 + G_0*(np.exp(packet.J_in/J_sat)-1))
                # if packet.J_in > 0 and J_out == 0: # J_in is too low, so we could approximate with small signal gain
                #     J_out = G_0 * packet.J_in
            # calculate output energy
            E_out = J_out * packet.A_in
            # calculate energy remain
            self.E_rem += (packet.E_in - E_out)

            # calculate B integral under construction
            average_intensity = (packet.E_in + E_out)/ (2 * packet.dfdt * packet.A_in)
            B_integral = 2 * pi * self.crystal.n2 * average_intensity * (self.crystal.length)/ (packet.wavelength)
            B_integral += packet.B_integral # add to the initial B integral in the signal

            # pack the result back to signal packet and add to output stack
            output_stack.append(SignalPacket(packet.wavelength, E_out, packet.A_in, packet.dfdt, packet.phase, B_integral))
        self.output_stack = np.array(output_stack)
        return self.output_stack

    # return output after amplification
    def get_output(self):
        return self.output_stack

    # return store object with store energy equal to remain energy
    def get_remain_energy(self):
        return Store(self.E_rem, self.store.A_sto)

class DataParser():
    def __init__(self, file_name):
        self.file_name = file_name
        self.data_frame = read_csv(file_name, sep = ',', index_col = 0)
        self.data_frame = self.data_frame.fillna('')

    def parse_unit_conversion(self, trial, parameter, si_unit):
        read = self.data_frame.loc[parameter, trial]
        unit = self.data_frame.loc[parameter, 'unit']
        return convert_unit(read, unit, si_unit)

    def parse_data(self, trial):
        # parse crystal data
        alpha = self.parse_unit_conversion(trial, 'crystal attenuation coefficient', '1/m')
        qd = self.parse_unit_conversion(trial, 'quantum defect', '')
        length = self.parse_unit_conversion(trial, 'crystal length', 'm')
        n2 = self.parse_unit_conversion(trial, 'nonlinear index', 'm^2/W')
        try:
            J_sat = self.parse_unit_conversion(trial, 'saturation fluence', 'J/m^2')
        except:
            J_sat = 0
        crystal = Crystal(alpha, qd, length, n2, J_sat)

        # parse pump data
        E_pump = self.parse_unit_conversion(trial, 'pump energy', 'J')
        A_pump = self.parse_unit_conversion(trial, 'pump area', 'm^2')
        pump = Pump(E_pump, A_pump)

        # parse initial signal data
        wavelength = self.parse_unit_conversion(trial, 'wavelength', 'm')
        pulse_duration = self.parse_unit_conversion(trial, 'pulse duration', 's')
        E_in = self.parse_unit_conversion(trial, 'initial input energy', 'J')
        A_in = self.parse_unit_conversion(trial, 'initial input area', 'm^2')
        initial_signal = SignalPacket(wavelength, pulse_duration, E_in, A_in, 4)

        return crystal, pump, initial_signal

class MultiPassAmplifier():
    # constructor
    def __init__(self, crystal, crystal_res = 1):
        # variables for simulation
        self.crystal = crystal
        self.crystal_res = crystal_res #1 is the same as not split the crystal
        self.loss_per_path = 0
        self.passes = []

    # delete last pass
    def remove_last_passes(self):
        self.passes.pop(-1)

    # delete all recorded passes
    def clear_passes(self):
        self.passes = []

    # return number of passes saved
    def get_passes_count(self):
        return len(self.passes)

    # set resolution for crystal splitting
    def set_crystal_res(self, res):
        if res >= 1: #if the given resolution is invalid, this is unchanged
            self.crystal_res = res
        return self.crystal_res

    # set percent loss per one path in crystal
    def set_loss_per_path(self, val):
        if val <= 1 and val >= 0: #if the given resolution is invalid, this is unchanged
            self.loss_per_path = val
        return self.loss_per_path

    # calculate attenuation to attenuation
    def absorbance_to_attenuation (self, length, absorbance, reflectance = 0):
        return -(1 / length) * log(1 - absorbance - reflectance)

    # use Beer - Lambert law to calculate transmittance
    def get_transmittance(self, attenuation, length):
        return np.exp(-attenuation * length)

    # calculate stored energy in crystal using Beer - Lambert law
    def calculate_stored_energy_from_pump(self, pump, dist = 'beer', front_pump = True, back_pump = True):
        store = []
        if not front_pump and not back_pump:
            raise Exception('No pump from neither front nor back side of the crystal')

        if dist == 'beer': # use Beer - Lambert law (longitudial pump)
            dl = self.crystal.length / self.crystal_res
            for i in range(self.crystal_res):
                #calculate effective pump energy
                f_pump = 0
                b_pump = 0
                if front_pump:
                    f_pump = pump.E_pump * self.get_transmittance(self.crystal.attenuation, i*dl) # front pump
                if back_pump:
                    b_pump = pump.E_pump * self.get_transmittance(self.crystal.attenuation, (self.crystal_res - i - 1)*dl) #back pump
                eff_E_pump = f_pump + b_pump
                #calculate absorbance
                abs = 1 - np.exp(-self.crystal.attenuation * dl) #assume reflectance = 0
                #calculate stored energy
                E_sto = eff_E_pump * self.crystal.qd * abs
                #calculate stored area
                A_sto = pump.A_pump
                #create a list of stored energy in crystal
                store.append(Store(E_sto, A_sto))
        elif dist == 'uniform' or dist == 'uni': # energy is uniformly distibuted (traverse pump)
            if front_pump and back_pump:
                eff_E_pump = 2 * pump.E_pump
            abs = 1 - np.exp(-self.crystal.attenuation * self.crystal.length) #assume reflectance = 0
            E_sto = eff_E_pump * self.crystal.qd * abs / self.crystal_res
            A_sto = pump.A_pump
            store = [Store(E_sto, A_sto) for _ in range(self.crystal_res)]
        else:
            raise Exception('Error pump energy distribution type')
        return np.array(store)

    # calculate single pass amplification
    def calculate_single_pass_amplification(self, store, signal):
        dl = self.crystal.length / self.crystal_res # length of each crystal split
        sto = store # current store energy in crystal
        in_sig = signal # current input signal pile
        input_train = [] # save input signal pile for each crystal
        output_train = [] # save output signal pile for each crystal
        amp_matrix = [] # save unit amplification for each crystal
        for i in range(self.crystal_res):
            # split crystal
            crys = self.crystal.copy()
            crys.length = dl
            amp_pile = []
            out_sig = [] # current output signal pile
            for sg in in_sig:
                # calculate output
                amp = UnitAmplifier(crys, sto[i], sg)
                output = amp.calculate_output()
                # save data for each crystal step
                amp_pile.append(amp)
                out_sig.append(output)
                # change stored energy to remain energy
                sto[i] = amp.get_remain_energy()
            out_sig = np.array(out_sig)
            # save data for each pulse step
            input_train.append(in_sig)
            output_train.append(out_sig)
            amp_matrix.append(amp_pile)
            # change output to input for next crystal split
            in_sig = out_sig
        input_train = np.array(input_train)
        output_train = np.array(output_train)
        amp_matrix = np.array(amp_matrix)
        self.passes.append((amp_matrix, input_train, output_train))
        return out_sig

    # calculate the final output signal coming out of the crystal
    def get_single_pass_output(self, path = -1):
        amp, input, output = self.passes[path]
        return output[-1]

    # calculate stored energy from remained energy after pass
    def calculate_stored_energy_from_remain(self):
        remain = []
        amp_matrix, *rest = self.passes[-1] #use data from last passes
        for amp_pile in amp_matrix:
            remain.append(amp_pile[-1].get_remain_energy()) #use remain energy after last time split passes
        return remain

    def multiply_signal_pile(self, signal_pile, multiplier):
        new_sig = []
        for signal_t in signal_pile:
            lst = []
            for signal_packet in signal_t:
                new_signal_pack = signal_packet.copy()
                new_signal_pack.E_in =  multiplier * new_signal_pack.E_in
                lst.append(new_signal_pack)
            new_sig.append(lst)
        return new_sig

    # calculate multi pass amplification
    def calculate_multi_pass_amplification(self, store, signal, run):
        sig = signal
        sto = store
        for i in range(run):
            sig = self.calculate_single_pass_amplification(sto, sig)
            sto = self.calculate_stored_energy_from_remain()
            sto.reverse()
            # account for loss per path
            if self.loss_per_path > 0 and self.loss_per_path <= 1:
                sig = self.multiply_signal_pile(sig, 1 - self.loss_per_path)
        return sig

    def get_total_energy_signal_pile(self, signal_pile):
        total_energy = 0
        for signal_t in signal_pile:
            for signal_packet in signal_t:
                total_energy += signal_packet.E_in
        return total_energy

    def calculate_single_pass_average_gain(self, path = -1):
        amp, input, output = self.passes[path]
        initial = input[0]
        final = output[-1]
        return self.get_total_energy_signal_pile(final) / self.get_total_energy_signal_pile(initial)

    def calculate_single_pass_average_B_integral(self, path = -1):
        amp, input, output = self.passes[path]
        initial = input[0]
        final = output[-1]
        n2 = amp[0][0].crystal.n2

    # feed current output back to amplification one time
    def feed_back_single_pass_amplification(self):
        amp_matrix, input_train, output_train = self.passes[path]
        last_signal_pile = output_train[-1]
        store = self.calculate_stored_energy_from_remain()
        return self.calculate_signal_pass_amplification(store, last_signal_pile)

    ### plot-related methods ###
    def get_crystal_length_list(self, crystal):
        dl = crystal.length / self.crystal_res
        return dl * np.arange(self.crystal_res)

    def plot_single_pass_store_energy(self, path = -1): #does not work
        amp, input, output = self.passes[path]
        rslt = []
        for amp_c in amp:
            rslt.append(amp_c[0].store.E_sto)
        rslt = np.array(rslt)
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_single_pass_remain_energy(self, path = -1):
        amp, input, output = self.passes[path]
        rslt = []
        for amp_c in amp:
            rslt.append(amp_c[-1].E_rem)
        rslt = np.array(rslt)
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_single_pass_input_energy(self, path = -1):
        amp, input, output = self.passes[path]
        rslt = np.array([self.get_total_energy_signal_pile(in_sig) for in_sig in input])
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_single_pass_output_energy(self, path = -1):
        amp, input, output = self.passes[path]
        rslt = np.array([self.get_total_energy_signal_pile(out_sig) for out_sig in output])
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_single_pass_gain(self, path = -1):
        amp, input, output = self.passes[path]
        ins = np.array([self.get_total_energy_signal_pile(in_sig) for in_sig in input])
        outs = np.array([self.get_total_energy_signal_pile(out_sig) for out_sig in output])
        rslt = outs/ins
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_single_pass_acc_gain(self, path = -1):
        amp, input, output = self.passes[path]
        initial = self.get_total_energy_signal_pile(input[0])
        outs = np.array([self.get_total_energy_signal_pile(out_sig) for out_sig in output])
        rslt = outs/initial
        crys = np.arange(self.crystal_res) * amp[0][0].crystal.length
        return crys, rslt

    def plot_multi_pass_parameter(self, single_pass_func):
        crys_tot = []
        rslt_tot = []
        for path in range(len(self.passes)):
            crys, rslt = single_pass_func(path)
            crys += path * self.crystal_res * (crys[1] - crys[0])
            crys_tot = np.concatenate((crys_tot, crys))
            rslt_tot = np.concatenate((rslt_tot, rslt))
        return crys_tot, rslt_tot

    def plot_multi_pass_store_energy(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_store_energy)

    def plot_multi_pass_remain_energy(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_remain_energy)

    def plot_multi_pass_input_energy(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_input_energy)

    def plot_multi_pass_output_energy(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_output_energy)

    def plot_multi_pass_gain(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_gain)

    def plot_multi_pass_acc_gain(self):
        return self.plot_multi_pass_parameter(self.plot_single_pass_acc_gain)
