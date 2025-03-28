import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
from fmPll import *
from fmSupportLib import *
from fmRRC import *
from fourierTransform import *
from scipy.signal import *

filter_taps = 101

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_decim = 10

audio_Fs = 48e3
mono_Fc = 16e3
mono_decim = 5

audio_Fs = 48e3

if __name__ == "__main__":
    rf_coeff = signal.firwin(filter_taps, rf_Fc/(rf_Fs/2), window = ('hann'))
    stereo_pilot_coeff = signal.firwin(filter_taps, [(18.5e3)/((rf_Fs/rf_decim)/2), (19.5e3/((rf_Fs/rf_decim)/2))], pass_zero="bandpass", window=('hann'))
    stereo_extract_coeff = signal.firwin(filter_taps, [(22e3)/((rf_Fs/rf_decim)/2), (54e3)/((rf_Fs/rf_decim)/2)], pass_zero="bandpass", window=('hann'))
    mono_coeff = signal.firwin(filter_taps, mono_Fc/((rf_Fs/rf_decim)/2), window=('hann'))
    stereo_lpf_coeff = signal.firwin(filter_taps, 38e3/((rf_Fs/rf_decim)/2), window=('hann'))
    # Print Bandpass Filter Coefficients
    print("Stereo Pilot BPF Coefficients:", stereo_pilot_coeff)
    print("Stereo Extract BPF Coefficients:", stereo_extract_coeff)

    # block processing variables
    block_size = 40e-3 * rf_Fs
    block_count = 0

    # Arrays to hold final audio data
    audio_data = np.array([])
    L = np.array([])
    R = np.array([])

    # set up subfigures for plotting
    subfig_height = np.array([0.8, 2, 1.6])
    plt.rc('figure', figsize=(7.5, 7.5))
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

    # read raw RF data from file
    in_fname = "data/iq_samples/stereo_l0_r9.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + str(len(raw_data)) + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
    bits = []
    print(len(iq_data))


    # State variables initialized to 0
    state_i = np.zeros(filter_taps-1)
    state_q = np.zeros(filter_taps-1)
    stereo_pilot_state = np.zeros(filter_taps-1)
    stereo_extract_state = np.zeros(filter_taps-1)
    pll_state = [0.0 ,0.0 ,1.0 ,0.0 ,1.0 ,0.0 ,0]
    stereo_state = np.zeros(filter_taps-1)
    mono_state = np.zeros(filter_taps-1)
    delay_state = np.zeros(int((filter_taps-1) / 2))
    # state_phase = 0.0
    state_phase = np.zeros(2)

    # Process the first 5 blocks
    # while block_count < 5 and (block_count + 1) * block_size * 2 <= len(iq_data):
    while (block_count + 1) * block_size * 2 <= len(iq_data):
        print(f"Block {block_count} - Data (first 10 samples): {iq_data[:10]}")

        # Extract a block of size block_size * 2 (to include both I and Q components)
        start_idx = block_count * block_size * 2
        end_idx = (block_count + 1) * block_size * 2
        iq_block = iq_data[int(start_idx):int(end_idx)]

        # Reshape into pairs of [I, Q]
        iq_block = np.reshape(iq_block, (-1, 2))

        # Split into I and Q components
        i_data = iq_block[:, 0]  # Extract I components
        q_data = iq_block[:, 1]  # Extract Q components

        print(f"Block {block_count} - I Data (first 10 samples): {i_data[:10]}")
        print(f"Block {block_count} - Q Data (first 10 samples): {q_data[:10]}")

        # LP filter the RF data
        i_rfFilter, state_i = signal.lfilter(rf_coeff, 1.0, i_data, zi=state_i)
        q_rfFilter, state_q = signal.lfilter(rf_coeff, 1.0, q_data, zi=state_q)

        # Downsample by 10
        i_data_ds = i_rfFilter[::rf_decim]
        q_data_ds = q_rfFilter[::rf_decim]

        print(f"Block {block_count} - Downsampled and LPF I Data (first 10 samples): {i_data_ds[:10]}")
        print(f"Block {block_count} - Downsampled and LPF Q Data (first 10 samples): {q_data_ds[:10]}")


        # Demodulate the RF data
        fm_demod, state_phase = fmDemodCustom(i_data_ds, q_data_ds, state_phase)

        # Print FM demodulated data
        print(f"Block {block_count} - FM Demodulated Data (first 10 samples): {fm_demod[:10]}")
        
        delay_fm_demod, delay_state = delayBlock(fm_demod, delay_state)

        # Extract pilot signal
        stereo_pilot, stereo_pilot_state = signal.lfilter(stereo_pilot_coeff, 1.0, fm_demod, zi=stereo_pilot_state)

        # Print pilot signal data
        print(f"Block {block_count} - Stereo Pilot Data (first 10 samples): {stereo_pilot[:10]}")

        # PLL to recover pilot tone
        stereo_pilotI, stereo_pilotQ,  pll_state = fmPll(stereo_pilot, 19e3, rf_Fs/rf_decim, pll_state, ncoScale=2.0)
        # plt.psd(stereo_pilotI, NFFT=1024, Fs=rf_Fs/10)
        # plt.show()

        # Extract stereo signal
        stereo_extract, stereo_extract_state = signal.lfilter(stereo_extract_coeff, 1.0, fm_demod, zi=stereo_extract_state)

        # Print stereo signal data
        print(f"Block {block_count} - Stereo Extract Data (first 10 samples): {stereo_extract[:10]}")

        # Mix stereo signal with pilot signal
        mixed_stereo = 2 * stereo_pilotI[:-1] * stereo_extract

        # Digital filtering and decimation
        stereo_predecim, stereo_state = signal.lfilter(stereo_lpf_coeff, 1.0, mixed_stereo, zi=stereo_state)
        stero = stereo_predecim[::mono_decim]

        # Extract mono
        mono_predecim, mono_state = signal.lfilter(mono_coeff, 1.0, delay_fm_demod, zi=mono_state)
        mono = mono_predecim[::mono_decim]

        # Stereo combiner
        L = np.concatenate((L, 0.5*(mono + stero)))
        R = np.concatenate((R, 0.5*(mono - stero)))

        print(f"Block {block_count} - L Data (first 10 samples): {L[-10:]}")
        print(f"Block {block_count} - R Data (first 10 samples): {R[-10:]}")

        block_count += 1

    # write left audio data to file
    out_fnameLeft = "data/fmStereoLeft.wav"
    wavfile.write(out_fnameLeft, int(audio_Fs), np.int16((L/2)*32767))
    print("Written audio samples to \"" + out_fnameLeft + "\" in signed 16-bit format")

    # write right audio data to file
    out_fnameRight = "data/fmStereoRight.wav"
    wavfile.write(out_fnameRight, int(audio_Fs), np.int16((R/2)*32767))
    print("Written audio samples to \"" + out_fnameRight + "\" in signed 16-bit format")

    # Combine left and right audio data in array [L, R]
    audio_data = np.column_stack((np.int16((L/2)*32767), np.int16((R/2)*32767)))

    # write audio data to file
    out_fname = "data/fmStereo.wav"
    wavfile.write(out_fname, int(audio_Fs), audio_data)
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
