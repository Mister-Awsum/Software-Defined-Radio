#include "Stereo.h"

// Function to implement a basic PLL for stereo carrier recovery
void fmPLL(const std::vector<dy4::real>& pllIn, 
           dy4::real freq, 
           dy4::real Fs, 
           dy4::real ncoScale, 
           dy4::real phaseAdjust, 
           dy4::real normBandwidth,
           std::vector<dy4::real>& pll_states,
           std::vector<dy4::real>& ncoOut) {
    // Scale factors for proportional/integrator terms (as in the Python version)
    dy4::real Cp = 2.666;
    dy4::real Ci = 3.555;

    // Gains for the loop filter
    dy4::real Kp = normBandwidth * Cp;
    dy4::real Ki = normBandwidth * normBandwidth * Ci;

    // Create output array for the NCO; size is one more than the input length
    ncoOut.clear();
    ncoOut.resize(pllIn.size()+1, 0.0);
    // Initialize internal state variables
    dy4::real integrator = pll_states[0];
    dy4::real phaseEst = pll_states[1];
    dy4::real feedbackI = pll_states[2];
    dy4::real feedbackQ = pll_states[3];
    ncoOut[0] = pll_states[4];
    dy4::real trigOffset = pll_states[5];
    dy4::real errorI =0.0;
    dy4::real errorQ = 0.0;
    dy4::real errorD = 0.0;
    dy4::real trigArg =0.0;

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
        trigOffset += 1;
        trigArg = 2 * M_PI * (freq / Fs) * trigOffset + phaseEst;
        feedbackI = cos(trigArg);
        feedbackQ = sin(trigArg);

        // Generate NCO output; note the addition of phaseAdjust (matching Python's phaseAdjust)
        ncoOut[k + 1] = cos(trigArg * ncoScale + phaseAdjust);
    }

    pll_states[0] = integrator;
    pll_states[1] = phaseEst;
    pll_states[2] = feedbackI;
    pll_states[3] = feedbackQ;
    pll_states[4] = ncoOut.back();
    pll_states[5] = trigOffset;
}

void delay_block(const std::vector<dy4::real> input_block, std::vector<dy4::real> &state_block, std::vector<dy4::real> &output_block) {    
    output_block.clear();

    for (size_t i = 0; i < state_block.size(); i++) {
        output_block.push_back(state_block[i]);
    }
    
    for (size_t i = 0; i < input_block.size() - state_block.size(); i++) {
        output_block.push_back(input_block[i]);
    }
    
    state_block.assign(input_block.end() - state_block.size(), input_block.end());
}
