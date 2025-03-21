/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "dy4.h"         // Defines dy4::real (e.g., typedef double)
#include "MonoBlock.h"   // Contains declarations for fmDemodArctan and fmDemodulate
#include "gtest/gtest.h"

namespace {

class FM_Demod_Fixture : public ::testing::Test {
public:
    const int N = 1024;  // Number of samples in our test signal
    const double EPSILON = 1e-6;  // Tolerance for floating-point comparisons

    // Test signals for I and Q channels
    std::vector<double> i_signal;
    std::vector<double> q_signal;
    
    // Outputs for demodulation
    std::vector<double> fm_demod_reference; // Computed via fmDemodArctan (reference)
    std::vector<double> fm_demod_custom;    // Computed via fmDemodulate (custom)

    // Phase state variables:
    double state_phase_ref = 0.0;
    std::vector<double> state_phase_custom;

    FM_Demod_Fixture() {
        i_signal.resize(N);
        q_signal.resize(N);
        fm_demod_reference.resize(N, 0.0);
        fm_demod_custom.resize(N, 0.0);
        state_phase_custom.resize(2, 0.0); // Custom implementation uses I/Q state vector
    }

    void generateFMSignal(std::vector<double>& signal, double min_val, double max_val) {
        for (int i = 0; i < N; i++) {
            signal[i] = min_val + (max_val - min_val) * i / (N - 1);
        }
    }

    void SetUp() override {
        // Generate test signals for I and Q
        generateFMSignal(i_signal, -1.0, 1.0);
        generateFMSignal(q_signal, -1.0, 1.0);

        // Ensure no value is too close to zero (to avoid division by zero)
        for (int i = 0; i < N; i++) {
            i_signal[i] = std::max(i_signal[i], 1e-6);
            q_signal[i] = std::max(q_signal[i], 1e-6);
        }

        // Compute reference FM demodulation using the arctan method.
        fmDemodArctan(i_signal, q_signal, state_phase_ref, fm_demod_reference);
    }

    void TearDown() override {}
};

TEST_F(FM_Demod_Fixture, FM_Demod_Arctan_NEAR) {
    std::vector<double> output(N, 0.0);
    double phase_state = state_phase_ref;
    fmDemodArctan(i_signal, q_signal, phase_state, output);

    ASSERT_EQ(fm_demod_reference.size(), output.size()) 
        << "Reference output and test output sizes differ.";

    for (int i = 0; i < (int)fm_demod_reference.size(); ++i) {
        EXPECT_NEAR(fm_demod_reference[i], output[i], EPSILON)
            << "Reference FM Demod mismatch at index " << i;
    }
}

TEST_F(FM_Demod_Fixture, FM_Demod_Custom_NEAR) {
    fmDemodulate(i_signal, q_signal, state_phase_custom, fm_demod_custom);

    ASSERT_EQ(fm_demod_reference.size(), fm_demod_custom.size()) 
        << "Reference and custom outputs have different sizes.";

    for (int i = 0; i < (int)fm_demod_reference.size(); ++i) {
        EXPECT_NEAR(fm_demod_reference[i], fm_demod_custom[i], EPSILON)
            << "FM Demod mismatch at index " << i;
    }

    EXPECT_NEAR(state_phase_custom[0], i_signal[N - 1], EPSILON)
        << "Custom state I value mismatch";
    EXPECT_NEAR(state_phase_custom[1], q_signal[N - 1], EPSILON)
        << "Custom state Q value mismatch";
}

} // namespace
