#include "RFFrontEnd.h"
#include "SyncQueue.h"

// Read data from stdin
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<dy4::real> &block_data) {
    // block_data.resize(num_samples, 0.0);

    std::vector<char> raw_data(num_samples);
    std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

    for(int k=0; k<(int)num_samples; k++) {
        block_data[k] = dy4::real(((unsigned char)raw_data[k] - 128) / 128.0);
    }
}

// Split I/Q samples from interleaved block data
void splitIQSamples(const std::vector<dy4::real> &block_data, std::vector<dy4::real> &i_samples, std::vector<dy4::real> &q_samples, unsigned int u_read_size) {
    for (int i = 0; i < u_read_size - 1; i += 2) {
        i_samples[i/2] = block_data[i];
        q_samples[i/2] = block_data[i+1];
    }
}

// RF Processing Thread
void RFProcessingThread(SyncQueue &syncQueue,
                        const std::vector<dy4::real> &rf_coeff,  int if_decim,
                        unsigned int u_read_size, int rf_taps) {
    std::vector<dy4::real> state_phase(2, 0.0); // State for FM demodulation
    std::vector<dy4::real> block_data(u_read_size);
    std::vector<dy4::real> state_i_lpf_100k(rf_taps - 1, 0.0);
    std::vector<dy4::real> state_q_lpf_100k(rf_taps - 1, 0.0);

    for (unsigned int block_id = 0; /*block_id<1*/; block_id++) {
        readStdinBlockData(u_read_size, block_id, block_data);

        if((std::cin.rdstate()) != 0) {
            std::cerr << "End of input stream reached" << std::endl;
            break;
        }

        // Split I/Q samples
        std::vector<dy4::real> i_samples(u_read_size / 2);
        std::vector<dy4::real> q_samples(u_read_size / 2);
        splitIQSamples(block_data, i_samples, q_samples, u_read_size);

        // RF Front-End Processing
        std::vector<dy4::real> i_ds(i_samples.size() / if_decim);
        std::vector<dy4::real> q_ds(q_samples.size() / if_decim);

        convolveAndDecimate(i_samples, rf_coeff, state_i_lpf_100k, if_decim, i_ds);
        convolveAndDecimate(q_samples, rf_coeff, state_q_lpf_100k, if_decim, q_ds);

        // FM Demodulation
        std::vector<dy4::real> fm_demod(i_ds.size());
        fmDemodulate(i_ds, q_ds, state_phase, fm_demod);

        // Push demodulated data to the output queue
        syncQueue.push({fm_demod, {}});
    }

    // Signal that RF processing is done
    syncQueue.setDone();
}

// FM DEMODULATION
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
    state_phase[0] = i_ds.back();
    state_phase[1] = q_ds.back();
}

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
