#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import cmath, math
import sys

def plotSpectrum(x, Fs, type = 'FFT'):

	n = len(x)             # length of the signal
	df = Fs / n            # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)
	XMag = abs(Xfreq) / n

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = XMag[0 : int(n / 2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs / 2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel = 'Frequency (Hz)', ylabel = 'Magnitude',
		   title = 'Frequency domain plot')
	# fig.savefig("freq.png")
	plt.show()

def plotTime(x, time):

	fig, ax = plt.subplots()
	ax.plot(time, x)
	ax.set(xlabel = 'Time (sec)', ylabel = 'Amplitude',
		   title = 'Time domain plot')
	# fig.savefig("time.png")
	plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

	dt = 1.0 / Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)      # time vector over interval

	# generate the sin signal
	x = amplitude * np.sin(2 * math.pi * frequency * time + phase)

	return time, x

def generateSquare(Fs, interval, frequency = 10.0, amplitude = 5.0, dutyCycle = 0.5):

	dt = 1.0 / Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)      # time vector over interval

	# generate the square signal
	x = amplitude * signal.square(2 * math.pi * frequency * time, duty = dutyCycle)

	return time, x

def fourierUnitTest(min = -1, max = 1, N = 32):

	# use a NumPy function to generate a random sequence of length "size"
	# of uniformly distributed numbers between min and max
	x = np.random.uniform(low = min, high = max, size = N)

	# the reconstructed signal (rx) must be the inverse Fourier transform of the frequency domain representation of x
	# below is a placeholder that helps us set up the unit test - eventually rx must be the real part of IDFT(DFT(x))
	rx = idft(dft(x))

	# check if all the values are close to each within the given tolerance
	if not np.allclose(x, rx, atol = 1e-4):
		print(f"Comparison between original signal and the inverse Fourier transform of its frequency domain representation fails")
		print(f"Original signal:", x)
		print(f"Reconstructed signal:", rx)
		print("Maximum difference:", np.max(np.abs(x - rx)))
	else:
		print(f"Unit test for Fourier/Inverse Fourier transform passed.")

def signalEnergyUnitTest(min = -10, max = 10, N = 1000):
	x = np.random.uniform(low = min, high = max, size = N)
	x_freq = dft(x)
	
	time_energy = signal_energy_time_domain(x)
	freq_energy = signal_energy_freq_domain(x_freq)
	
	if not np.isclose(time_energy, freq_energy, atol = 1e-4):
		print(f"Comparison between time domain energy and frequency domain energy fails")
		print(f"Time domain energy:", time_energy)
		print(f"Frequency domain energy:", freq_energy)
	else:
		print(f"Unit test for Parseval's theorem passed.")

def dft(x):
	x_freq = []

	for m in range(len(x)):
		x_freq.append(0)
		for k in range(len(x)):
			x_freq[m] += x[k] * cmath.exp(-2*math.pi*1j*((k*m)/len(x)))
	return x_freq

def idft(x):
	x_time = []
	for m in range(len(x)):
		x_time.append(0)
		for k in range(len(x)):
			x_time[m] += (x[k] * cmath.exp(2*math.pi*1j*((k*m)/len(x)))) / len(x)
	return x_time

def signal_energy_time_domain(x_time):
	signalEnergy = 0
	for k in range(len(x_time)):
		signalEnergy += abs(x_time[k]) ** 2
	return signalEnergy

def signal_energy_freq_domain(x_freq):
	signalEnergy = 0
	for k in range(len(x_freq)):
		signalEnergy += (abs(x_freq[k]) ** 2) / len(x_freq)
	return signalEnergy
	
	
def generateTritoneSignal(Fs, interval):
	sin1_time, sin1_x = generateSin(Fs, interval, 3.0, 5.0, 0.0)
	sin2_time, sin2_x = generateSin(Fs, interval, 7.0, 2.0, 20.0)
	sin3_time, sin3_x = generateSin(Fs, interval, 14.0, 10.0, 30.0)
	
	return sin1_time, (sin1_x + sin2_x + sin3_x)

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\til3: in-lab 3')
	print('\tth:  take-home')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0          # sampling rate
	interval = 1.0      # set up to one full second

	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for the Fourier transform')

		# generate the user-defined sin function
		time, x = generateSin(Fs, interval)
		# plot the signal in time domain
		plotTime(x, time)
		# plot the signal in frequency domain
		plotSpectrum(x, Fs, type = 'FFT')

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the Fourier transform')

		# compute the spectrum with your own DFT
		# you can use cmath.exp() for complex exponentials
		# plotSpectrum(x, Fs, type = 'your DFT name')
		"""time, x = generateSin(Fs, interval)
		x_dft = dft(x)
		x_idft = idft(x_dft)
		plotSpectrum(x_dft, Fs, type = 'FFT')
		plotTime(x_idft, time)"""
		# confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
		# the self check should be done within the fourierUnitTest function
		fourierUnitTest()
		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the Fourier transform')

		# you can use np.random.randn() to generate a random parameter
		# we can overwrite the default values
		"""frequency = 8.0                      # frequency of the signal
		amplitude = 3.0                      # amplitude of the signal
		phase = 1.0                          # phase of the signal
		time, x = generateSin(Fs, interval, frequency, amplitude, phase)
		x_freq = dft(x)
		
		freq_energy = signal_energy_freq_domain(x_freq)
		time_energy = signal_energy_time_domain(x)"""
		# You should also numerically check if the signal energy
		# in time and frequency domains is identical
		#
		# create your own Unit Test to validate Parseval's theorem
		# (based on the same principle as for the Fourier transform)
		signalEnergyUnitTest()
		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'il3'):

		print('In-lab experiment 3 for the Fourier transform')

		# generate randomized multi-tone signals
		time, x = generateTritoneSignal(Fs, interval)
		# plot them in both time and frequency domain
		x_freq = dft(x)
		plotTime(x, time)
		plotSpectrum(x_freq, Fs, type = 'FFT')

		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for the Fourier transform')

		# for specific details check the lab document
		time, x = generateSquare(Fs, interval)
		plotTime(x, time)
		plotSpectrum(x, Fs, type = 'FFT')
	else:

		cli_error_msg()

	plt.show()
