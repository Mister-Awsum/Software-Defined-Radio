#include "Stereo.h"

// Function to implement a basic PLL for stereo carrier recovery
void fmPLL(const std::vector<dy4::real>& pllIn, 
           dy4::real freq, 
           dy4::real Fs, 
           dy4::real ncoScale, 
           dy4::real phaseAdjust, 
           dy4::real normBandwidth, 
           std::vector<dy4::real>& output) {
    // Scale factors for proportional/integrator terms (as in the Python version)
    dy4::real Cp = 2.666;
    dy4::real Ci = 3.555;

    // Gains for the loop filter
    dy4::real Kp = normBandwidth * Cp;
    dy4::real Ki = pow(normBandwidth, 2) * Ci;

    // Create output array for the NCO; size is one more than the input length
    std::vector<dy4::real> ncoOut(pllIn.size() + 1, 0.0);

    // Initialize internal state variables
    dy4::real integrator = 0.0;
    dy4::real phaseEst = 0.0;
    dy4::real feedbackI = 1.0;
    dy4::real feedbackQ = 0.0;
    dy4::real trigOffset = 0.0;
    dy4::real errorI;
    dy4::real errorQ;
    dy4::real errorD;
    dy4::real trigArg;

    ncoOut[0] = 1.0;

    // Process each sample in pllIn
    for (size_t k = 0; k < pllIn.size(); k++) {
        // Phase detector: compute error using the complex conjugate of the feedback signal
        errorI = pllIn[k] * feedbackI;
        errorQ = pllIn[k] * (-feedbackQ);

        // Four-quadrant arctangent discriminator for phase error detection
        errorD = atan2(errorQ, errorI);

        // Loop filter: update the integrator and phase estimate
        integrator += Ki * errorD;
        phaseEst += Kp * errorD + integrator;

        // Internal oscillator: update trigOffset and compute trigArg
        trigOffset += 1.0;
        trigArg = 2 * M_PI * (freq / Fs) * trigOffset + phaseEst;
        feedbackI = cos(trigArg);
        feedbackQ = sin(trigArg);

        // Generate NCO output; note the addition of phaseAdjust (matching Python's phaseAdjust)
        ncoOut[k + 1] = cos(trigArg * ncoScale + phaseAdjust);
    }

    // Return the NCO output
    output = ncoOut;
}

void s_processing(vector<dy4::real> s_carrier, vector<dy4::real> s_channel, vector<short int> &audio_data) {
    // Step 1: Downconvert the stereo channel by mixing with the recovered carrier.
    size_t N = std::min(s_carrier.size(), s_channel.size());

    vector<dy4::real> stereo_diff(N);
    for (size_t i = 0; i < N; i++) {
        stereo_diff[i] = s_channel[i] * s_carrier[i];
    }

    // Step 2: Convert mono audio data (L+R) from 16-bit integers to floating point.
    size_t M = audio_data.size();
    if (M == 0) {
        // cerr << "Warning: audio_data is empty before processing!" << endl;
        return;
    }

    vector<dy4::real> mono_real(M);
    for (size_t i = 0; i < M; i++) {
        mono_real[i] = static_cast<dy4::real>(audio_data[i]) / 16384.0; 
    }

    // Step 3: Ensure stereo_diff and mono_real are the same size
    size_t minSize = std::min(mono_real.size(), stereo_diff.size());
    if (minSize == 0) {
        // cerr << "Warning: stereo_diff or mono_real has zero size!" << endl;
        return;
    }

    // Step 4: Compute Left and Right Channels
    vector<dy4::real> left(minSize, 0.0);
    vector<dy4::real> right(minSize, 0.0);
    for (size_t i = 0; i < minSize; i++) {
        left[i] = 0.5 * (mono_real[i] + stereo_diff[i]);
        right[i] = 0.5 * (mono_real[i] - stereo_diff[i]);
    }

    // Step 5: Convert left and right channels to interleaved stereo output.
    vector<short int> stereo_audio(2 * minSize, 0);
    for (size_t i = 0; i < minSize; i++) {
        stereo_audio[2 * i] = static_cast<short int>(left[i] * 16384);      // Left channel
        stereo_audio[2 * i + 1] = static_cast<short int>(right[i] * 16384); // Right channel
    }

    // Step 6: Assign to audio_data only if valid
    if (!stereo_audio.empty()) {
        audio_data = stereo_audio;
    } else {
        // cerr << "Warning: stereo_audio is empty after processing!" << endl;
    }
}


// (Am / 2) * (cos(2 * M_PI * (fc + fm) * t) + cos(2 * M_PI * (fc - fm) * t));
// (Am / 4) * (2 * cos(2 * M_PI * fm * t) + cos(2 * M_PI * (fc + fm) * t) + cos(2 * M_PI * (fc - fm) * t));