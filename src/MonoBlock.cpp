#include "MonoBlock.h"

// Function to upsample / downsample a signal
void upsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &upsampled) {
    if (factor <= 1) {
        upsampled = signal;
        return;
    }
    upsampled.resize(signal.size() * factor, 0.0);  // Allocate output
    
    for (size_t i = 0; i < signal.size(); i++) {
        upsampled[i * factor] = signal[i];
    }
}

void downsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &downsampled) {
    if (factor <= 1) {
        downsampled = signal;
        return;
    }
    downsampled.resize(signal.size() / factor, 0.0);  // Allocate output

    for (size_t i = 0; i < downsampled.size(); i++) {
        downsampled[i] = signal[i * factor];
    }
}

// FM DEMODULATION
void fmDemodArctan(const std::vector<dy4::real> &i_ds, const std::vector<dy4::real> &q_ds, dy4::real &state_phase, std::vector<dy4::real> &fm_demod) {
    fm_demod.resize(i_ds.size());

    for (size_t n = 0; n < i_ds.size(); n++) {
        // Compute phase using atan2
        dy4::real phase_current = atan2(q_ds[n], i_ds[n]);

        // Compute phase difference
        dy4::real phase_diff = phase_current - state_phase;

        // Phase unwrapping (keep within -π to π)
        if (phase_diff > M_PI) {
            phase_diff -= 2 * M_PI;
        } else if (phase_diff < -M_PI) {
            phase_diff += 2 * M_PI;
        }

        fm_demod[n] = phase_diff;

        // Update previous phase for continuity
        state_phase = phase_current;
    }
}

void fmDemodulate(const std::vector<dy4::real> &i_ds, const std::vector<dy4::real> &q_ds, std::vector<dy4::real> &state_phase, std::vector<dy4::real> &fm_demod) {  
    for (size_t n = 0; n < i_ds.size(); n++) {
 
        // For the previous sample, if we are at n==0, use the saved state
        dy4::real i_prev = (n > 0) ? i_ds[n - 1] : state_phase[0];
        dy4::real q_prev = (n > 0) ? q_ds[n - 1] : state_phase[1];


        // Compute the differences 
        dy4::real i_delta = i_ds[n] - i_prev;   // dI
        dy4::real q_delta = q_ds[n] - q_prev;   // dQ
 
        // Instead of computing a difference of arctan's, we use the differential method:
 
        //   fm_demod[n] = (i[n] * (q_next - q_prev) - (i_next - i_prev) * q[n]) / (2*(i[n]^2 + q[n]^2))
        dy4::real prod_1 = i_ds[n] * q_delta;
        dy4::real prod_2 = i_delta * q_ds[n];
        dy4::real denom = pow(i_ds[n], 2) + pow(q_ds[n], 2);
 
        fm_demod[n] = (denom != 0) ? ((prod_1 - prod_2) / denom) : 0.0;
    }
 
    // Save the last I and Q samples for state continuity in the next block
    state_phase[0] = i_ds[i_ds.size() - 1];
    state_phase[1] = q_ds[q_ds.size() - 1];
}

void generate_fm_values(std::vector<dy4::real> &signal, dy4::real min_val, dy4::real max_val) {
    for (auto& val : signal) {
        val = min_val + (max_val - min_val) * (rand() / (RAND_MAX + 1.0));
    }
}

// Read data from stdin
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<dy4::real> &block_data) {
    // block_data.resize(num_samples, 0.0);

    std::vector<char> raw_data(num_samples);
    std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

    for(int k=0; k<(int)num_samples; k++) {
        block_data[k] = dy4::real(((unsigned char)raw_data[k] - 128) / 128.0);
    }
}

// Generate FIR filter coefficients
// void firwin(int num_taps, dy4::real cutoff, dy4::real fs, vector<dy4::real> &coeffs) {
//     coeffs.resize(num_taps, 0.0);
//     dy4::real norm_cutoff = cutoff / (fs / 2);
//     for (int i = 0; i < num_taps; i++) {
//         dy4::real n = i - (num_taps - 1) / 2.0;
//         if (n == 0.0) {
//             coeffs[i] = norm_cutoff;
//         } else {
//             coeffs[i] = norm_cutoff * sin(M_PI * norm_cutoff * n) / (M_PI * norm_cutoff * n);
//         }
//         coeffs[i] *= 0.5 - 0.5 * cos(2.0 * M_PI * i / (num_taps - 1)); // Hann window
//     }
//     dy4::real sum = accumulate(coeffs.begin(), coeffs.end(), 0.0);
//     for (auto& coeff : coeffs) coeff /= sum; // Normalize filter
// }

// Function to process I/Q data in blocks
// void processIQData(const vector<dy4::real> &iq_data, vector<dy4::real> &audio_data) {
//     int block_size = 1024 * rf_decim * audio_decim * 2;
//     int block_count = 0;

//     // Filter coefficients
//     vector<real> rf_coeff, audio_coeff;
//     firwin(rf_taps, rf_Fc, rf_Fs, push_backrf_coeff);
//     firwin(audio_taps, 15e3, audio_Fs, audio_coeff); // Example: 15 kHz cutoff

//     // State variables
//     vector<real> state_i_lpf_100k(rf_taps - 1, 0.0);
//     vector<real> state_q_lpf_100k(rf_taps - 1, 0.0);
//     vector<real> state_audio(audio_taps - 1, 0.0);
//     real state_phase = 0.0;

//     while (static_cast<size_t>((block_count + 1) * block_size) < iq_data.size()) {
//         cerr << "Processing block " << block_count << endl;

//         // Extract I/Q samples
//         vector<real> i_samples, q_samples;
//         for (int i = block_count * block_size; i < (block_count + 1) * block_size; i += 2) {
//             i_samples.push_back(iq_data[i]);
//             q_samples.push_back(iq_data[i + 1]);
//         }

//         // Apply low-pass filter
//         vector<real> i_filt, q_filt;
//         applyFIR(i_samples, rf_coeff, state_i_lpf_100k, i_filt);
//         applyFIR(q_samples, rf_coeff, state_q_lpf_100k, q_filt);

//         // Downsample
//         downsample(i_filt, rf_decim);
//         downsample(q_filt, rf_decim);

//         // FM Demodulation
//         vector<real> fm_demod;
//         fmDemodulate(i_filt, q_filt, state_phase, fm_demod);

//         // Apply audio filtering
//         vector<real> audio_filt;
//         applyFIR(fm_demod, audio_coeff, state_audio, audio_filt);

//         // Downsample audio
//         vector<real> audio_block;
//         downsample(audio_filt, audio_decim);

//         // Append to audio buffer
//         audio_data.insert(audio_data.end(), audio_filt.begin(), audio_filt.end());

//         block_count++;
//     }

//     cerr << "Finished processing all the blocks from the recorded I/Q samples" << endl;
// }

// Function to write audio data to a WAV file
// void writeWAV(const std::string &filename, const std::vector<dy4::real> &audio_data, const dy4::real sample_rate) {
//     std::ofstream file(filename, std::ios::binary);
//     if (!file) {
//         std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
//         return;
//     }

//     // Convert floating-point samples to 16-bit PCM
//     std::vector<int16_t> pcm_data(audio_data.size());
//     for (size_t i = 0; i < audio_data.size(); i++) {
//         pcm_data[i] = static_cast<int16_t>(std::max(-1.0, std::min(1.0, audio_data[i])) * 32767);
//     }

//     // WAV Header
//     uint32_t data_chunk_size = pcm_data.size() * sizeof(int16_t);
//     uint32_t chunk_size = 36 + data_chunk_size;

//     file.write("RIFF", 4);
//     file.write(reinterpret_cast<const char *>(&chunk_size), 4);
//     file.write("WAVE", 4);

//     // Format chunk
//     file.write("fmt ", 4);
//     uint32_t fmt_chunk_size = 16;
//     uint16_t audio_format = 1; // PCM
//     uint16_t num_channels = 1; // Mono
//     uint32_t sample_rate_int = static_cast<uint32_t>(sample_rate);
//     uint16_t bits_per_sample = 16;
//     uint16_t block_align = num_channels * (bits_per_sample / 8);
//     uint32_t byte_rate = sample_rate_int * block_align;

//     file.write(reinterpret_cast<const char *>(&fmt_chunk_size), 4);
//     file.write(reinterpret_cast<const char *>(&audio_format), 2);
//     file.write(reinterpret_cast<const char *>(&num_channels), 2);
//     file.write(reinterpret_cast<const char *>(&sample_rate_int), 4);
//     file.write(reinterpret_cast<const char *>(&byte_rate), 4);
//     file.write(reinterpret_cast<const char *>(&block_align), 2);
//     file.write(reinterpret_cast<const char *>(&bits_per_sample), 2);

//     // Data chunk
//     file.write("data", 4);
//     file.write(reinterpret_cast<const char *>(&data_chunk_size), 4);
//     file.write(reinterpret_cast<const char *>(pcm_data.data()), data_chunk_size);

//     file.close();
// }

// void normalizeVector(std::vector<real> &vec, real new_min, real new_max) {
// 	real max_val = *std::max_element(vec.begin(), vec.end());
// 	real min_val = *std::min_element(vec.begin(), vec.end());

// 	real original_range = max_val - min_val;
// 	real new_range = new_max - new_min;
//  push_back
// 	for (int i = 0; i < vec.size(); i++) {
// 		vec[i] = new_min + (vec[i] - min_val) * new_range / original_range;
// 	}

// }
