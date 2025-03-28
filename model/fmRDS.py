import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from fmPll import fmPll
from fmSupportLib import *
from fmRRC import *
from fourierTransform import *
from scipy.signal import *
from fmDataProcessing import process_rds_from_constellation

# === Configuration Parameters ===
NUM_TAPS = 101
RF_SAMPLE_RATE = 2.4e6
RF_CENTER_FREQ = 100e3
DECIMATION = 10

# Mode 0 RDS parameters
UPSAMPLE_FACTOR = 133
DOWNSAMPLE_FACTOR = 480
# Effective sample rate for RRC filter (~66.5 kHz)
EFFECTIVE_RRC_RATE = (RF_SAMPLE_RATE / DECIMATION) * (UPSAMPLE_FACTOR / DOWNSAMPLE_FACTOR)

# === Helper Functions (Rewritten) ===

def shift_signal_samples(sig_array, prev_hist):
    shifted = np.concatenate((prev_hist, sig_array[:-len(prev_hist)]))
    new_hist = sig_array[-len(prev_hist):]
    return shifted, new_hist

def calc_windowed_fir(f_low, f_high, s_rate, num_coeff):
    center = (num_coeff - 1) / 2.0
    mid_freq = (f_low + f_high) / 2.0
    coeffs = np.zeros(num_coeff)
    norm_sum = 0.0
    for k in range(num_coeff):
        n = k - center
        if n == 0:
            coeffs[k] = (2 * f_high / s_rate) - (2 * f_low / s_rate)
        else:
            coeffs[k] = (np.sin(2 * np.pi * f_high * n / s_rate) / (np.pi * n)) - (np.sin(2 * np.pi * f_low * n / s_rate) / (np.pi * n))
        coeffs[k] *= (0.5 - 0.5 * np.cos(2 * np.pi * k / (num_coeff - 1)))
        norm_sum += coeffs[k] * np.cos(2 * np.pi * n * mid_freq / s_rate)
    return coeffs / norm_sum

def fir_resample(ds_factor, up_factor, x, blk_len, fir_coeff, state):
    out_len = ((len(x) * up_factor) + ds_factor - 1) // ds_factor
    y = [0.0] * out_len
    for i in range(out_len):
        up_idx = i * ds_factor
        phase = up_idx % up_factor
        sample_val = 0.0
        for j in range((len(fir_coeff) + phase) // up_factor):
            coeff_idx = j * up_factor + phase
            if coeff_idx >= len(fir_coeff):
                continue
            x_idx = up_idx // up_factor - j
            if 0 <= x_idx < len(x):
                sample_val += fir_coeff[coeff_idx] * x[x_idx]
            elif x_idx < 0 and abs(x_idx) <= len(state):
                sample_val += fir_coeff[coeff_idx] * state[len(state) + x_idx]
        y[i] = sample_val
    new_state = np.zeros(len(state))
    history_needed = (len(fir_coeff) + up_factor - 1) // up_factor - 1
    for j in range(min(history_needed, len(x))):
        new_state[len(state) - j - 1] = x[len(x) - j - 1]
    return y, new_state

def recover_timing(i_vals, q_vals, sym_samples, gain=0.1, err_prev=0.0, phase_prev=0.0):
    i_vals = np.array(i_vals)
    q_vals = np.array(q_vals)
    rec_i, rec_q = [], []
    phase = phase_prev
    err = err_prev
    idx = 0
    last_sample = None
    while True:
        pos = int(idx + phase)
        rec_i.append(i_vals[pos])
        rec_q.append(q_vals[pos])
        if last_sample is not None:
            mid = int(idx + phase - sym_samples/2)
            if 0 <= mid < len(i_vals):
                mid_val = i_vals[mid]
                avg_val = (i_vals[pos] + last_sample) / 2.0
                err = (avg_val - mid_val) * gain
            phase += err
            while phase >= sym_samples:
                phase -= sym_samples
            while phase < 0:
                phase += sym_samples
        last_sample = i_vals[pos]
        if pos + sym_samples >= len(i_vals):
            break
        idx += sym_samples
    rem_i = i_vals[idx:] if idx < len(i_vals) else np.array([])
    rem_q = q_vals[idx:] if idx < len(q_vals) else np.array([])
    return rec_i, rec_q, err, phase, rem_i, rem_q

def manchester_decode(samples, last_val, init_phase, err_count):
    decoded = []
    if init_phase == 1:
        if last_val > 0 and samples[0] < 0:
            decoded.append(1)
        elif last_val < 0 and samples[0] > 0:
            decoded.append(0)
        else:
            decoded.append(-1)
            err_count += 1
    phase = init_phase
    while phase + 1 < len(samples):
        if samples[phase] > 0 and samples[phase+1] < 0:
            decoded.append(1)
        elif samples[phase] < 0 and samples[phase+1] > 0:
            decoded.append(0)
        else:
            decoded.append(-1)
            err_count += 1
        phase += 2
        if phase+1 >= len(samples):
            break
    print("Manchester block length:", len(decoded), "Initial phase:", init_phase)
    init_phase = 0 if init_phase == 1 else 1
    return decoded, samples[-1], init_phase, err_count

def diff_decode(bits, init_state=0):
    out = [init_state ^ bits[0]]
    for i in range(1, len(bits)):
        out.append(bits[i-1] ^ bits[i])
    return out, bits[-1]

# === End of Helper Functions ===

# === Filter Designs ===
rf_filter = signal.firwin(NUM_TAPS, RF_CENTER_FREQ/(RF_SAMPLE_RATE/2), window='hann')
rds_extract_filter = signal.firwin(NUM_TAPS, [54e3/((RF_SAMPLE_RATE/DECIMATION)/2), 60e3/((RF_SAMPLE_RATE/DECIMATION)/2)],
                                   pass_zero="bandpass", window='hann')
rds_pilot_filter = signal.firwin(NUM_TAPS, [113.5e3/((RF_SAMPLE_RATE/DECIMATION)/2), 114.5e3/((RF_SAMPLE_RATE/DECIMATION)/2)],
                                 pass_zero="bandpass", window='hann')
rds_lowpass = signal.firwin(NUM_TAPS, 3e3/((RF_SAMPLE_RATE/DECIMATION)/2), window='hann')
rcos_filter = impulseResponseRootRaisedCosine(EFFECTIVE_RRC_RATE, NUM_TAPS)

temp_fs = RF_SAMPLE_RATE/DECIMATION
resampFir = signal.firwin(NUM_TAPS * UPSAMPLE_FACTOR, min((UPSAMPLE_FACTOR/DOWNSAMPLE_FACTOR)*temp_fs, temp_fs/2)/((temp_fs*UPSAMPLE_FACTOR)/2), window='hann')
resampFir = [val * UPSAMPLE_FACTOR for val in resampFir]

# === Initialize State Variables ===
state_I = np.zeros(NUM_TAPS-1)
state_Q = np.zeros(NUM_TAPS-1)
extract_state = np.zeros(NUM_TAPS-1)
pilot_state = np.zeros(NUM_TAPS-1)
pll_state = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0]
pre_resamp_I = np.zeros(NUM_TAPS-1)
pre_resamp_Q = np.zeros(NUM_TAPS-1)
rrc_state_I = np.zeros(NUM_TAPS-1)
rrc_state_Q = np.zeros(NUM_TAPS-1)
state_phase = 0.0
delay_memory = np.zeros(int((NUM_TAPS-1)/2))
man_phase = 0
last_man_sample = 0
err_count = 0
diff_state = 0

# IMPORTANT: Define both resampler states
resample_state_I = np.zeros(NUM_TAPS * UPSAMPLE_FACTOR - 1)
resample_state_Q = np.zeros(NUM_TAPS * UPSAMPLE_FACTOR - 1)

# === Block Processing Setup ===
blk_size = int(480000)  # Reduced block size (e.g., ~480,000 samples ~200ms of RF data)
blk_count = 0
constellation_I = []
constellation_Q = []
diff_bit_stream = []

# === Read Raw IQ Data ===
in_filename = "data/iq_samples/samples3.raw"
raw_data = np.fromfile(in_filename, dtype='uint8')
print("Read raw RF data from \"" + str(len(raw_data)) + "\" in unsigned 8-bit format")
iq_data = (np.float32(raw_data) - 128.0)/128.0
print("Converted IQ data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
print("Total IQ samples:", len(iq_data))

# === Main Processing Loop ===
while (blk_count+1)*blk_size < len(iq_data):
    print("Processing block", blk_count)
    I_data = iq_data[blk_count*blk_size:(blk_count+1)*blk_size:2]
    Q_data = iq_data[blk_count*blk_size+1:(blk_count+1)*blk_size:2]
    
    filt_I, state_I = signal.lfilter(rf_filter, 1.0, I_data, zi=state_I)
    filt_Q, state_Q = signal.lfilter(rf_filter, 1.0, Q_data, zi=state_Q)
    
    ds_I = filt_I[::DECIMATION]
    ds_Q = filt_Q[::DECIMATION]
    
    fm_out, state_phase = fmDemodArctan(ds_I, ds_Q, state_phase)
    
    rds_out, extract_state = signal.lfilter(rds_extract_filter, 1.0, fm_out, zi=extract_state)
    rds_sq = rds_out * rds_out
    pilot_out, pilot_state = signal.lfilter(rds_pilot_filter, 1.0, rds_sq, zi=pilot_state)
    
    # Correct fmPll call: pass entire pll_state as fourth parameter (7 parameters total)
    nco_out, nco_Q, pll_state = fmPll(pilot_out, 114e3, (RF_SAMPLE_RATE/DECIMATION), pll_state, 0.5, 0, 0.002)
    
    delayed_rds, delay_memory = shift_signal_samples(rds_out, delay_memory)
    
    pilot_complex = nco_out + 1j * nco_Q
    pilot_complex = pilot_complex[:-1]
    baseband = delayed_rds * np.conj(pilot_complex)
    mixed_I = np.real(baseband)
    mixed_Q = np.imag(baseband)
    
    preI, pre_resamp_I = signal.lfilter(rds_lowpass, 1.0, mixed_I, zi=pre_resamp_I)
    preQ, pre_resamp_Q = signal.lfilter(rds_lowpass, 1.0, mixed_Q, zi=pre_resamp_Q)
    
    resamp_I, resample_state_I = fir_resample(DOWNSAMPLE_FACTOR, UPSAMPLE_FACTOR, preI, blk_size, resampFir, resample_state_I)
    resamp_Q, resample_state_Q = fir_resample(DOWNSAMPLE_FACTOR, UPSAMPLE_FACTOR, preQ, blk_size, resampFir, resample_state_Q)
    
    rrc_I, rrc_state_I = signal.lfilter(rcos_filter, 1.0, resamp_I, zi=rrc_state_I)
    rrc_Q, rrc_state_Q = signal.lfilter(rcos_filter, 1.0, resamp_Q, zi=rrc_state_Q)
    
    blk_count += 1
    
    # Symbol extraction: fixed window (28 samples per symbol)
    sps = 28
    num_syms = 2375
    half_win = 8
    for sym in range(num_syms):
        center = sym * sps
        start = center - half_win
        end = center + half_win
        if start < 0 or end >= len(rrc_I):
            continue
        seg = np.abs(rrc_I[start:end])
        peak = np.argmax(seg)
        peak_idx = start + peak
        constellation_I.append(rrc_I[peak_idx])
        constellation_Q.append(rrc_Q[peak_idx])
    
    # After processing enough blocks, perform Manchester and differential decoding
    if blk_count > 20:
        # Try decoding with both possible initial phases
        man_bits0, last0, phase0, err0 = manchester_decode(constellation_I, last_man_sample, 0, 0)
        man_bits1, last1, phase1, err1 = manchester_decode(constellation_I, last_man_sample, 1, 0)
        # Choose the decoding with fewer errors
        if err0 <= err1:
            man_bits, last_man_sample, man_phase, err_count = man_bits0, last0, phase0, err0
        else:
            man_bits, last_man_sample, man_phase, err_count = man_bits1, last1, phase1, err1

        diff_bits_out, diff_state = diff_decode(man_bits, diff_state)
        diff_bit_stream.extend([int(b) for b in diff_bits_out])
    
    print(f"Block {blk_count}: RRC samples = {len(rrc_I)}, symbols extracted = {len(constellation_I)}")
    
# Normalize final constellation for display
constellation_I = np.array(constellation_I)
constellation_Q = np.array(constellation_Q)
constellation_I -= np.mean(constellation_I)
constellation_Q -= np.mean(constellation_Q)
norm_val = np.max(np.sqrt(constellation_I**2 + constellation_Q**2))
constellation_I /= norm_val
constellation_Q /= norm_val

plt.scatter(constellation_I, constellation_Q, s=1)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Final Constellation Diagram")
plt.show()

# --- Clean the bitstream --- 
# Since errors produce negative values, clip them to produce a binary sequence.
clean_bits = [0 if b < 0 else 1 for b in diff_bit_stream]

# Process the recovered differential bitstream into RDS messages
result = process_rds_from_constellation(clean_bits)
