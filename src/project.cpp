// Authors: Brandon English, Deyontae Patterson, Marco, Victoria Black
// Date: 03/27/2025
// Description: Main file for the project. This file reads in raw IQ data from stdin, processes it, and outputs audio data to stdout.

#include "dy4.h"
#include "iofunc.h"
#include "logfunc.h"
#include "MonoBlock.h"
#include "Stereo.h"

using namespace std;
using namespace dy4;

// THE MAIN FUNCTION
int main(int argc, char *argv[]) {
    const dy4::real rf_Fc = 100e3;  // RF cutoff frequency
    const dy4::real if_Fc = 16e3;  // IF cutoff frequency
    const int rf_taps = 101;        // Number of taps for LPF
    const int if_taps = 101;
    const int audio_taps = 101;     // Number of taps for audio filter

    dy4::real pll_freq;
    dy4::real pll_phase;
    dy4::real normBandwidth;
    dy4::real nco_Scalar;

    if (string(argv[1]) == "s") {
        pll_freq = 19e3;
        pll_phase = 0.0;
        normBandwidth = 0.01;
        nco_Scalar = 2.0;
    } else if (string(argv[1]) == "r") {
        pll_freq = 0.0;
        pll_phase = 0.0;
        normBandwidth = 0.01;
        nco_Scalar = 0.5;
    }

    // Variables dependent on Mode
    dy4::real audio_Fs;
    dy4::real if_Fs;
    dy4::real rf_Fs;
    int BLOCK_SIZE;
    int U;

    if (string(argv[2]) == "0") {
        audio_Fs = 48e3;
        if_Fs = 240e3;
        rf_Fs = 2400e3;
        BLOCK_SIZE = 96000;
        U = 1;
    } else if (string(argv[2]) == "1") {
        audio_Fs = 32e3;
        if_Fs = 128e3;
        rf_Fs = 1152e3;
        BLOCK_SIZE = 46080;
        U = 1;
    } else if (string(argv[2]) == "2") {
        audio_Fs = 44.1e3;
        if_Fs = 240e3;
        rf_Fs = 2400e3;
        BLOCK_SIZE = 96000;
        U = 147;
    } else {
        audio_Fs = 44.1e3;
        if_Fs = 288e3;
        rf_Fs = 2304e3;
        BLOCK_SIZE = 92160;
        U = 49;
    }

    int if_decim = rf_Fs / if_Fs;
    int D = if_Fs * U / audio_Fs;

    cerr << "IF Decim = " << if_decim << ", Upsampling factor U = " << U << ", Downsampling Factor D = " << D << endl;

    // Filter coefficients
    vector<dy4::real> rf_coeff(rf_taps);
    vector<dy4::real> if_coeff(if_taps*U);
    vector<dy4::real> extract_coeff(if_Fs);
    vector<dy4::real> pilot_coeff(if_Fs);
    vector<dy4::real> rds_coeff(if_Fs);
    vector<dy4::real> demod_coeff(if_Fs);
    vector<dy4::real> rds_demod_coeff(if_Fs);
    vector<dy4::real> carrier_coeff(if_Fs);
    vector<dy4::real> state_coeff(if_Fs);

    generateLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
    generateLPF(if_Fs*U, if_Fc, if_taps*U, if_coeff);

    for(int n=0; n<if_coeff.size(); n++) {
        if_coeff[n] *= U;
    }

    generateBPF(pilot_coeff, 18.5e3, 19.5e3, if_Fs, if_taps);
    generateBPF(extract_coeff, 22e3, 54e3, if_Fs, if_taps);

    generateBPF(rds_coeff, 54e3, 60e3, if_Fs, if_taps);
    generateBPF(carrier_coeff, 113.5e3, 114.5e3, if_Fs, if_taps);

    // State variables
    vector<dy4::real> state_i_lpf_100k(rf_taps - 1, 0.0);
    vector<dy4::real> state_q_lpf_100k(rf_taps - 1, 0.0);
    vector<dy4::real> state_audio(if_taps*U - 1, 0.0);
    vector<dy4::real> state_extract(if_taps - 1, 0.0);
    vector<dy4::real> state_recovery(if_taps - 1, 0.0);
    vector<dy4::real> state_rds(if_taps - 1, 0.0);
    vector<dy4::real> state_phase(2, 0.0);

    // Initialize Vector S izes
    dy4::real read_size = BLOCK_SIZE * 2;
    dy4::real iq_size = read_size / 2;
    dy4::real demod_size = iq_size / if_decim;
    dy4::real demod_up_size = demod_size * U;
    dy4::real audio_size = demod_up_size / D;

    cerr << "Read Size = " << read_size << endl;
    cerr << "IQ Size = " << iq_size << endl;
    cerr << "Demod Size = " << demod_size << endl;
    cerr << "Audio Size = " << audio_size << endl;

    unsigned int u_read_size = (unsigned int)read_size;

    vector<dy4::real> block_data(u_read_size);

    vector<dy4::real> i_samples((unsigned int)iq_size);
    vector<dy4::real> q_samples((unsigned int)iq_size);

    vector<dy4::real> i_ds((unsigned int)demod_size);
    vector<dy4::real> q_ds((unsigned int)demod_size);
    vector<dy4::real> fm_demod((unsigned int)demod_size);
    vector<dy4::real> demod_up((unsigned int)demod_up_size);

    vector<dy4::real> processed_data((unsigned int)audio_size);
    vector<short int> audio_data((unsigned int)audio_size);

    vector<dy4::real> path_bpf;
    vector<dy4::real> rds_recovered;
    vector<dy4::real> rds_carrier;
    vector<dy4::real> rds_filt;
    vector<dy4::real> state_carrier;

    // Other

    // MAIN PROGRAM LOOP
    for(unsigned int block_id = 0; /*block_id<1*/; block_id++) {
        readStdinBlockData(u_read_size, block_id, block_data);

        if((std::cin.rdstate()) != 0) {
            std::cerr << "End of input stream reached" << std::endl;
            break;
        }

        for (int i=0; i<u_read_size-1; i+=2) {
            i_samples[i/2] = block_data[i];
            q_samples[i/2] = block_data[i+1];
        }

        // FRONT-END PROCESSING
        convolveAndDecimate(i_samples, rf_coeff, state_i_lpf_100k, if_decim, i_ds);
        convolveAndDecimate(q_samples, rf_coeff, state_q_lpf_100k, if_decim, q_ds);

        fmDemodulate(i_ds, q_ds, state_phase, fm_demod);

        if(block_id == 10) {writeBinData("../data/fm_demod.bin", fm_demod);}

        // MONO PATH
        if (string(argv[1]) == "m" || string(argv[1]) == "s") {
            
            // upsample(fm_demod, U, demod_up);

            // convolveAndDecimate(demod_up, if_coeff, state_audio, D, processed_data);

            ConvolveFast(fm_demod, if_coeff, state_audio, U, D, processed_data);

            printRealVector(state_audio);
 
            for(unsigned int k=0; k<processed_data.size(); k++) {
                if (std::isnan(processed_data[k])) {
                    cerr << "Audio at index " << k << " set to 0" << endl;
                    audio_data[k] = 0;  
                } 

                else audio_data[k] = static_cast<short int>(processed_data[k] * 16384);
            }

            // STEREO PATH
            if(string(argv[1]) == "s") {
                vector<short int> stereo_audio = audio_data;

                // Stereo Channel Extraction
                vector<dy4::real> s_channel(fm_demod.size());
                convolveAndDecimate(fm_demod, extract_coeff, state_extract, 1, s_channel);

                // Stereo Carrier Recovery
                vector<dy4::real> stereo_filt(fm_demod.size());
                convolveAndDecimate(fm_demod, pilot_coeff, state_recovery, 1, stereo_filt);

                vector<dy4::real> s_carrier(fm_demod.size());
                fmPLL(stereo_filt, pll_freq, audio_Fs, nco_Scalar, pll_phase, normBandwidth, s_carrier);
                
                // Stereo Processing
                s_processing(s_carrier, s_channel, stereo_audio);

                // [???] Downsample final audio
                // for (size_t i = 0; i < audio_data.size(); i++) {
                //     audio_data[i] = stereo_audio[i * 2];
                // }
            }

            fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
        }

        // RADIO MODE
        else if (string(argv[1]) == "r") {
            // RDS CHANNEL EXTRACTION
            vector<dy4::real> rds_channel(fm_demod.size());
            convolveBlock(fm_demod, rds_coeff, state_rds, rds_channel);

            // RDS CARRIER RECOVERY
            vector<dy4::real> rds_path = rds_channel;

            for(int i=0; i < rds_path.size(); i++) {
                rds_path[i] = pow(rds_path[i], 2);
            }

            convolveBlock(rds_path, carrier_coeff, state_carrier, path_bpf);

            fmPLL(path_bpf, pll_phase, audio_Fs, nco_Scalar, pll_phase, normBandwidth, rds_recovered);

            /* [TO-DO] Implement **All Filters block*/

            // RDS DEMODULATION
            size_t N = std::min(rds_carrier.size(), rds_channel.size());

            vector<dy4::real> rds_diff(N);
            for (size_t i = 0; i < N; i++) {
                rds_diff[i] = rds_channel[i] * rds_carrier[i];
            }

            convolveBlock(rds_diff, demod_coeff, rds_demod_coeff, rds_filt);

            // RDS DATA PROCESSING

        }

    }

    const std::string in_fname = "../data/fm_demod.bin";
	std::vector<dy4::real> bin_data;
	readBinData(in_fname, bin_data);

	// Generate an index vector to be used by logVector on the X axis
	std::vector<dy4::real> vector_index;
	genIndexVector(vector_index, bin_data.size());
	logVector("demod_time", vector_index, bin_data);

	// Take a slice of data with a limited number of samples for the Fourier transform
	std::vector<dy4::real> slice_data = std::vector<dy4::real>(bin_data.begin(), bin_data.begin() + NFFT);

	std::vector<std::complex<dy4::real>> Xf;
	DFT(slice_data, Xf);

	std::vector<dy4::real> Xmag;
	computeVectorMagnitude(Xf, Xmag);

	// Log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // Log only positive freq

	vector_index.clear();
	estimatePSD(slice_data, if_Fs, vector_index, Xmag);
	logVector("demod_psd", vector_index, Xmag);

	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);

    // Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n

    return 0;
}

// TESTS
// int main(int argc, char *argv[]) {
//     const dy4::real rf_Fc = 100e3;  // RF cutoff frequency
//     const int rf_taps = 101;        // Number of taps for LPF
//     const int if_taps = 101;
//     const int audio_taps = 101;     // Number of taps for audio filter

//     // Variables dependent on Mode
//     dy4::real audio_Fs;
//     dy4::real if_Fs;
//     dy4::real rf_Fs;
//     dy4::real block_factor;

//     if (string(argv[2]) == "0") {
//         audio_Fs = 48e3;
//         if_Fs = 240e3;
//         rf_Fs = 2.4e6;
//         block_factor = 1;
//     } else if (string(argv[2]) == "1") {
//         audio_Fs = 32e3;
//         if_Fs = 128e3;
//         rf_Fs = 1.152e6;
//         block_factor = 0.48;
//     } else if (string(argv[2]) == "2") {
//         audio_Fs = 44.1e3;
//         if_Fs = 240e3;
//         rf_Fs = 2.4e6;
//         block_factor = 1;
//     } else {
//         audio_Fs = 44.1e3;
//         if_Fs = 288e3;
//         rf_Fs = 2.304e6;
//         block_factor = 0.96;
//     }

//     const dy4::real rf_decim = rf_Fs / if_Fs;
//     const dy4::real audio_decim = if_Fs / audio_Fs;
//     const int BLOCK_SIZE = 96000 * block_factor;

//     // Filter coefficients
//     vector<dy4::real> rf_coeff, audio_coeff, extract_coeff, pilot_coeff, dc_block_coeff;

//     generateLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
//     generateLPF(audio_Fs, 16e3, audio_taps, audio_coeff);

//     // generateBPF(pilot_coeff, 18.5e3, 19.5e3, audio_Fs, audio_taps);
//     // generateBPF(extract_coeff, 22e3, 54e3, audio_Fs, audio_taps);

//     // State variables
//     vector<dy4::real> state_i_lpf_100k(rf_taps - 1, 0.0);
//     vector<dy4::real> state_q_lpf_100k(rf_taps - 1, 0.0);
//     // vector<dy4::real> state_audio(audio_taps - 1, 0.0);
//     // vector<dy4::real> state_extract(audio_taps - 1, 0.0);
//     // vector<dy4::real> state_recovery(audio_taps - 1, 0.0);
//     vector<dy4::real> state_phase(2, 0.0);

//     std::vector<dy4::real> block_data(BLOCK_SIZE, 0.0);

//     int block_id = 0;
//     readStdinBlockData(BLOCK_SIZE, block_id, block_data);

//     if((std::cin.rdstate()) != 0) {
//         std::cerr << "End of input stream reached" << std::endl;
//         exit(1);
//     }

//     // Extract I/Q samples
//     vector<dy4::real> i_samples(BLOCK_SIZE/2, 0.0);
//     vector<dy4::real> q_samples(BLOCK_SIZE/2, 0.0);
//     for (int i=0; i<BLOCK_SIZE-1; i+=2) {
//         i_samples[i/2] = block_data[i];
//         q_samples[(i-1)/2] = block_data[i+1];
//     }     
        
//     vector<dy4::real> i_filt, q_filt;
//     convolveBlock(i_samples, rf_coeff, state_i_lpf_100k, i_filt);
//     convolveBlock(q_samples, rf_coeff, state_q_lpf_100k, q_filt);

//     downsample(i_filt, rf_decim);
//     downsample(q_filt, rf_decim);

//     vector<dy4::real> fm_demod(i_filt.size(), 0.0);
//     // vector<dy4::real> fm_demod(i_samples.size(), 0.0);
//     fmDemodulate(i_filt, q_filt, state_phase, fm_demod);
//     // fmDemodulate(i_samples, q_samples, state_phase, fm_demod);

//     writeBinData("../data/fm_demod.bin", fm_demod);

//     const std::string in_fname = "../data/fm_demod.bin";
// 	std::vector<dy4::real> bin_data;
// 	readBinData(in_fname, bin_data);

// 	// Generate an index vector to be used by logVector on the X axis
// 	std::vector<dy4::real> vector_index;
// 	genIndexVector(vector_index, bin_data.size());
// 	logVector("demod_time", vector_index, bin_data);

// 	// Take a slice of data with a limited number of samples for the Fourier transform
// 	std::vector<dy4::real> slice_data = std::vector<dy4::real>(bin_data.begin(), bin_data.begin() + NFFT);

// 	std::vector<std::complex<dy4::real>> Xf;
// 	DFT(slice_data, Xf);

// 	std::vector<dy4::real> Xmag;
// 	computeVectorMagnitude(Xf, Xmag);

// 	// Log the frequency magnitude vector
// 	vector_index.clear();
// 	genIndexVector(vector_index, Xmag.size());
// 	logVector("demod_freq", vector_index, Xmag); // Log only positive freq

// 	vector_index.clear();
// 	estimatePSD(slice_data, if_Fs, vector_index, Xmag);
// 	logVector("demod_psd", vector_index, Xmag);

// 	const std::string out_fname = "../data/outdata.bin";
// 	writeBinData(out_fname, bin_data);

//     // Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n

//     return 0;
// }