#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

def custom_firwin(num_taps, cutoff, sample_rate, window='hann'):
    nyquist = sample_rate / 2
    normalized_cutoff = cutoff / nyquist

    # Create a sinc filter
    n = np.arange(0, num_taps) - (num_taps - 1) / 2
    h = np.sinc(normalized_cutoff * n)

    # Apply the window function
    if window == 'hann':
        w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(num_taps) / (num_taps - 1))
        h *= w

    # Normalize the filter coefficients
    h /= np.sum(h)
    return h.astype(np.float32)

def custom_convolve(filter_coefficients, feedback_coefficient, input_signal):
    # Ensure feedback_coefficients[0] is normalized to 1
    if feedback_coefficient != 1.0:
        filter_coefficients = filter_coefficients / feedback_coefficient

    # Ensure proper padding
    filter_len = len(filter_coefficients)
    padded_input = np.pad(input_signal, (filter_len - 1, 0), mode='constant')

    output_signal = np.zeros(len(input_signal), dtype=np.float64)

    for n in range(len(input_signal)):
        output_signal[n] = np.sum(filter_coefficients * padded_input[n:n + filter_len])

    return output_signal

def blockConvolveFIR(input_block, h, state, blockSize):
    # Prepend previous state to input block
    input_with_state = np.concatenate((state, input_block))

    # Perform convolution
    output_with_state = custom_convolve(h, 1.0, input_with_state)

    # Extract the current output block
    output_block = output_with_state[len(state):len(state) + blockSize]

    # Update state for next block (store last (len(h) - 1) samples)
    new_state = input_with_state[-(len(h) - 1):]
    output_block = np.array(output_block, dtype=np.float32)

    return output_block, new_state

def custom_fmDemod(i, q, state=None):
    if state is None:
        state = [0.0, 0.0]

    i_state, q_state = state[0], state[1]
    
    demodulated = np.zeros(len(i))  # Initialize output array
    
    for n in range(len(i)):
        # Compute the differentials
        i_prev = i[n - 1] if n > 0 else i_state
        q_prev = q[n - 1] if n > 0 else q_state
        
        i_next = i[n + 1] if n + 1 < len(i) else i[n]
        q_next = q[n + 1] if n + 1 < len(q) else q[n]
        
        i_delta = (i_next - i_prev)
        q_delta = (q_next - q_prev)
        
        # Compute the product of the differentials
        prod_1 = i[n] * q_delta
        prod_2 = i_delta * q[n]
        
        # Compute denominator
        denom = 2 * (i[n] ** 2 + q[n] ** 2)
        
        # final division
        demodulated[n] = (prod_1 - prod_2) / denom if denom != 0 else 0.0
        
    # Update state
    state = [i[-1], q[-1]]
    
    return demodulated, state

rf_Fs = 1152e3
rf_Fc = 100e3
rf_taps = 101
rf_decim = 9

if_Fs = 128e3

audio_Fs = 32e3
audio_decim = 4
# add other settings for audio, like filter taps, ...
audio_taps = 101
audio_Fc = 16e3
# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs take-home (il_vs_th = 1)
il_vs_th = 1

if il_vs_th == 0:
    audio_coeff = signal.firwin(audio_taps, audio_Fc / (if_Fs / 2), window='hann')
else:
    audio_coeff = custom_firwin(audio_taps, audio_Fc, if_Fs, window='hann')

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/iq_samples/samples0.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	'''
	# IQ data is normalized between -1 and +1 in 64-bit double format
	iq_data = (np.float64(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 64-bit double format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
	'''

	# coefficients for the front-end low-pass filter
	if(il_vs_th == 0):
		rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))
	else:
		rf_coeff = custom_firwin(rf_taps, cutoff=rf_Fc, sample_rate=rf_Fs, window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(audio_taps, audio_Fc / (if_Fs / 2), window='hann')
	else:
		# to be updated by you for the take-home exercise
		# with your own code for impulse response generation
		audio_coeff = custom_firwin(audio_taps, cutoff=audio_Fc, sample_rate=if_Fs, window='hann')

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 46080

	if (il_vs_th == 0):
		# in-lab implementation
		block_count = 0
  
 		# states needed for continuity in block processing
		state_i_lpf_100k = np.zeros(rf_taps - 1)
		state_q_lpf_100k = np.zeros(rf_taps - 1)
		state_phase = 0
		state_audio = np.zeros(audio_taps - 1)
		# add state as needed for the mono channel filter

		# audio buffer that stores all the audio blocks
		audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

		# if the number of samples in the last block is less than the block size
		# it is fine to ignore the last few samples from the raw IQ file
		while (block_count + 1) * block_size < len(iq_data):

			# if you wish to have shorter runtimes while troubleshooting
			# you can control the above loop exit condition as you see fit
			#print('Processing block ' + str(block_count))

			# filter to extract the FM channel (I samples are even, Q samples are odd)
			i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size:(block_count + 1) * block_size:2], zi=state_i_lpf_100k)
			q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2], zi=state_q_lpf_100k)

			# downsample the I/Q data from the FM channel
			i_ds = i_filt[::rf_decim]
			q_ds = q_filt[::rf_decim]

			# FM demodulator
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

			# extract the mono audio data through filtering
			# downsample audio data
			# audio_block = ... change as needed
			audio_filt, state_audio = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=state_audio)
			audio_block = audio_filt[::audio_decim]
			audio_data = np.concatenate((audio_data, audio_block))

			# concatenate the most recently processed audio_block
			# to the previous blocks stored already in audio_data
			#
			#

			# to save runtime, select the range of blocks to log data
			# this includes both saving binary files and plotting PSD
			if block_count >= 10 and block_count < 12:

				# plot PSD of selected block after FM demodulation
				# (for easier visualization purposes we divide Fs by 1e3 to imply the kHz units on the x-axis)
				# (this scales the y axis of the PSD, but not the relative strength of different frequencies)
				ax0.clear()
				fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], \
						'Demodulated FM (block ' + str(block_count) + ')')
				# output binary file name (where samples are written from Python)
				fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
				# create binary file where each sample is a 32-bit float
				fm_demod.astype('float32').tofile(fm_demod_fname)

				'''
				# create binary file where each sample is a 64-bit double
				fm_demod.astype('float64').tofile(fm_demod_fname)
				'''

				ax1.clear()
				fmPlotPSD(ax1, audio_filt, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Extracted Mono')

				ax2.clear()
				fmPlotPSD(ax2, audio_data, audio_Fs / 1e3, subfig_height[2], 'Downsampled Mono Audio')
	
				# save figure to file
				fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

			block_count += 1

	else:
		block_count = 0
  
		# Take-home exercise implementation
		state_i_lpf_100k = np.zeros(rf_taps - 1)
		state_q_lpf_100k = np.zeros(rf_taps - 1)
  
		state_phase = None
  
		state_audio = np.zeros(audio_taps - 1)

		audio_data = np.array([])

		while (block_count + 1) * block_size < len(iq_data):
			# Filtering I and Q channels
			#i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size:(block_count + 1) * block_size:2], zi=state_i_lpf_100k)
			#q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2], zi=state_q_lpf_100k)

			i_filt, state_i_lpf_100k = blockConvolveFIR(h=rf_coeff, input_block=iq_data[block_count * block_size:(block_count + 1) * block_size:2], state=state_i_lpf_100k, blockSize=block_size)
			q_filt, state_q_lpf_100k = blockConvolveFIR(h=rf_coeff, input_block=iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2], state=state_q_lpf_100k, blockSize=block_size)
   
			# Downsample I/Q data
			i_ds = i_filt[::rf_decim]
			q_ds = q_filt[::rf_decim]

			# FM demodulation (unchanged)
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, state_phase = custom_fmDemod(i_ds, q_ds, state_phase)

			audio_filt, state_audio = blockConvolveFIR(input_block=fm_demod, h=audio_coeff, state=state_audio, blockSize=block_size)

			# Downsample audio data
			audio_block = audio_filt[::audio_decim]
			audio_data = np.concatenate((audio_data, audio_block))

			# Save or plot if in the specified block range
			if block_count >= 10 and block_count < 12:
				ax0.clear()
				fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], 'Demodulated FM (block ' + str(block_count) + ')')

				fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
				fm_demod.astype('float32').tofile(fm_demod_fname)

				ax1.clear()
				fmPlotPSD(ax1, audio_filt, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Extracted Mono')

				ax2.clear()
				fmPlotPSD(ax2, audio_data, audio_Fs / 1e3, subfig_height[2], 'Downsampled Mono Audio')

				fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

			block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data / 2) * 32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	plt.show()
