#include "AudioProcessing.h"

void AudioProcessingThread(SyncQueue &syncQueue, const std::vector<dy4::real> &if_coeff, 
                           std::vector<dy4::real> &state_audio, int U, int D, 
                           const std::vector<dy4::real> &extract_coeff, std::vector<dy4::real> &state_extract, std::vector<dy4::real> &stereo_process_state,
                           const std::vector<dy4::real> &pilot_coeff, std::vector<dy4::real> &state_recovery, std::vector<dy4::real> &state_block, std::vector<dy4::real> &output_block,
                           dy4::real pll_freq, dy4::real pll_phase, dy4::real if_Fs, dy4::real normBandwidth, dy4::real nco_Scalar,
                           std::vector<dy4::real> &pll_states, const std::string &mode) {
    if (mode != "m" && mode != "s") {
        // Do nothing if mode is not "m" or "s"
        return;
    }

    while (true) {
        std::pair<std::vector<dy4::real>, std::vector<dy4::real>> rf_data;
        if (!syncQueue.pop(rf_data)) {
            break; // Exit if the queue is empty and done flag is set
        }

        const auto &fm_demod = rf_data.first;

        if (mode == "s") {
            // STEREO PROCESSING
            // Step 1: Stereo Carrier Recovery
            std::vector<dy4::real> audio_data(fm_demod.size() * U / D);
            std::vector<dy4::real> processed_data(fm_demod.size() * U / D);
            vector<dy4::real> stereo_filt(fm_demod.size());
            ConvolveFast(fm_demod, pilot_coeff, state_recovery, 1, 1, stereo_filt);
            //convolveAndDecimate(fm_demod, pilot_coeff, state_recovery, 1, stereo_filt);
            //convolveBlock(fm_demod, pilot_coeff, state_recovery, stereo_filt);
            vector<dy4::real> s_carrier(fm_demod.size());
            fmPLL(stereo_filt, pll_freq, if_Fs, nco_Scalar, pll_phase, normBandwidth, pll_states, s_carrier);

            // Step 2: Stereo Channel Extraction
            vector<dy4::real> s_channel(fm_demod.size());
            // convolveAndDecimate(fm_demod, extract_coeff, state_extract, 1, s_channel);
            ConvolveFast(fm_demod, extract_coeff, state_extract, 1, 1, s_channel);
            //convolveBlock(fm_demod, extract_coeff, state_extract, s_channel);
            
            // Step 3: Stereo mixing
            size_t N = std::min(s_carrier.size(), s_channel.size());
            vector<dy4::real> stereo_diff(N);
            for (size_t i = 0; i < N; i++) {
                stereo_diff[i] = 2 * s_channel[i] * s_carrier[i];
            }
            vector<dy4::real> mono_real(audio_data.size());
            size_t minSize = std::min(mono_real.size(), stereo_diff.size());

            // Step 4: Filter and decimate the stereo signal
            vector<dy4::real> stereo(minSize);
            ConvolveFast(stereo_diff, if_coeff, stereo_process_state, U, D, stereo);
            // convolveAndDecimate(stereo_diff, rf_coeff, stereo_process_state, D, stereo);

            // Step 5: Allpass filter and filter mono data
            delay_block(fm_demod, state_block, output_block);
            // convolveAndDecimate(output_block, if_coeff, state_audio, D, processed_data);
            ConvolveFast(output_block, if_coeff, state_audio, U, D, processed_data);
            // Step 6: Compute Left and Right Channels
            vector<dy4::real> left(minSize, 0.0);
            vector<dy4::real> right(minSize, 0.0);
            for (size_t i = 0; i < minSize; i++) {
                left[i] = 0.5 * (processed_data[i] + stereo[i]);
                right[i] = 0.5 * (processed_data[i] - stereo[i]);
            }

            // Step 6: Convert left and right channels to interleaved stereo output.
            vector<short int> stereo_audio(2 * minSize, 0);
            for (size_t i = 0; i < minSize; i++) {
                stereo_audio[2 * i] = static_cast<short int>(left[i] * 16384);      // Left channel
                stereo_audio[2 * i + 1] = static_cast<short int>(right[i] * 16384); // Right channel
            }

            fwrite(&stereo_audio[0], sizeof(short int), stereo_audio.size(), stdout);
        } else if (mode == "m") {
            // MONO PROCESSING
            std::vector<dy4::real> processed_data(fm_demod.size() * U / D);
            std::vector<short int> audio_data(fm_demod.size() * U / D);
            ConvolveFast(fm_demod, if_coeff, state_audio, U, D, processed_data);

            for (size_t i = 0; i < processed_data.size(); i++) {
                audio_data[i] = static_cast<short int>(processed_data[i] * 16384);
            }

            fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
        }
    }
}
