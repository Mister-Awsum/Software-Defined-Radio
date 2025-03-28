// Authors: Brandon English, Deyontae Patterson, Marco Sobol, Victoria Black
// Date: 03/27/2025
// Description: Main file for the project. This file reads in raw IQ data from stdin, processes it, and outputs audio data to stdout.

#include "dy4.h"
#include "iofunc.h"
#include "logfunc.h"
#include "MonoBlock.h"
#include "Stereo.h"
#include "RFFrontEnd.h"
#include <queue>
#include <mutex>
#include "SyncQueue.h"
#include <thread>
#include "AudioProcessing.h"

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
    bool isStereo;

    if (string(argv[1]) == "s") {
        pll_freq = 19e3;
        pll_phase = 0.0;
        normBandwidth = 0.01;
        nco_Scalar = 2.0;
        isStereo = true;
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
        if_Fs = 192e3;
        if_Fs = 192e3;
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
    vector<dy4::real> extract_coeff(if_taps);
    vector<dy4::real> pilot_coeff(if_taps);
    vector<dy4::real> demod_coeff(if_Fs);
    vector<dy4::real> carrier_coeff(if_Fs);
    vector<dy4::real> state_coeff(if_Fs);
    
    generateLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
    generateLPF(if_Fs*U, if_Fc, if_taps*U, if_coeff);

    // Upscale the intermediate frequency coefficients by a factor of U, explained in project document
    for(int n=0; n<if_coeff.size(); n++) {
        if_coeff[n] *= U;
    }

    generateBPF(pilot_coeff, 18.5e3 , 19.5e3, if_Fs, if_taps);
    generateBPF(extract_coeff, 22e3, 54e3, if_Fs, if_taps);

    // State variables
    vector<dy4::real> state_i_lpf_100k(rf_taps - 1, 0.0);
    vector<dy4::real> state_q_lpf_100k(rf_taps - 1, 0.0);
    vector<dy4::real> stereo_process_state(rf_taps - 1, 0.0);
    vector<dy4::real> state_audio(if_taps*U - 1, 0.0);
    vector<dy4::real> state_extract(if_taps - 1, 0.0);
    vector<dy4::real> state_recovery(if_taps - 1, 0.0);
    vector<dy4::real> state_phase(2, 0.0);
    vector<dy4::real> pll_states = {0, 0, 1, 0, 1, 0};

    // Initialize Vector Sizes
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

    vector<dy4::real> i_samples((unsigned int)iq_size);
    vector<dy4::real> q_samples((unsigned int)iq_size);

    vector<dy4::real> i_ds((unsigned int)demod_size);
    vector<dy4::real> q_ds((unsigned int)demod_size);
    vector<dy4::real> fm_demod((unsigned int)demod_size);
    vector<dy4::real> demod_up((unsigned int)demod_up_size);

    vector<dy4::real> output_block(fm_demod.size());
    vector<dy4::real> state_block((if_taps - 1) / 2, 0.0);

    vector<dy4::real> processed_data((unsigned int)audio_size);
    vector<short int> audio_data((unsigned int)audio_size);
    
    vector<dy4::real> mono_real;
    size_t minSize;
    vector<dy4::real> stereo;
    
    // Variables for multi-threading with a queue
    SyncQueue syncQueue; // Queue for RF-processed data

    
    // Start the RF processing thread for the current block
    std::thread rf_thread(RFProcessingThread, std::ref(syncQueue),
                            std::ref(rf_coeff), if_decim, u_read_size, rf_taps);

    // Start the audio processing thread
    std::thread audio_thread(AudioProcessingThread, std::ref(syncQueue), std::ref(if_coeff), 
                           std::ref(state_audio), U, D, 
                           std::ref(extract_coeff), std::ref(state_extract), std::ref(stereo_process_state),
                           std::ref(pilot_coeff), std::ref(state_recovery), std::ref(state_block), std::ref(output_block),
                           pll_freq, pll_phase, if_Fs, normBandwidth, nco_Scalar,
                           std::ref(pll_states), string(argv[1]));

    // Wait for the RF thread to finish
    rf_thread.join();

    // Wait for the audio thread to finish
    audio_thread.join();

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
