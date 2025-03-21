// filter_test.cpp
/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include "dy4.h"       // Defines dy4::real (e.g. float or double)
#include "filter.h"    // Declaration of generateLPF() and generateBPF()
#include "gtest/gtest.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace {

class Filter_Fixture : public ::testing::Test {
public:
    // For LPF testing
    const dy4::real Fs = 2.4e6;   // Example RF sampling rate
    const dy4::real Fc = 100e3;   // LPF cutoff frequency
    const int lpf_taps = 101;     // Number of taps for LPF
    std::vector<dy4::real> lpf;

    // For BPF testing
    const dy4::real fs = 48000;   // Example sample rate for audio (or IF)
    const dy4::real f_low = 22000;  // Lower cutoff frequency for BPF (22 kHz)
    const dy4::real f_high = 54000; // Higher cutoff frequency for BPF (54 kHz)
    const int bpf_taps = 101;       // Number of taps for BPF
    std::vector<dy4::real> bpf;

    Filter_Fixture() {}

    void SetUp() override {
        // Generate LPF and BPF coefficients using your functions
        generateLPF(Fs, Fc, lpf_taps, lpf);
        generateBPF(bpf, f_low, f_high, fs, bpf_taps);
    }

    void TearDown() override {}
};

// Test that the LPF is symmetric (i.e., h[k] == h[N-1-k])
TEST_F(Filter_Fixture, LPF_Symmetry) {
    for (int k = 0; k < lpf_taps; k++) {
        EXPECT_NEAR(lpf[k], lpf[lpf_taps - 1 - k], 1e-6)
            << "LPF symmetry mismatch at index " << k;
    }
}

// Test that the center tap of the LPF equals the normalized cutoff value
// (i.e. NormCutOff = Fc/(Fs/2) as used in generateLPF)
TEST_F(Filter_Fixture, LPF_CenterTap) {
    int mid = (lpf_taps - 1) / 2;
    dy4::real expected_center = Fc / (Fs / 2);
    EXPECT_NEAR(lpf[mid], expected_center, 1e-6)
        << "LPF center tap value incorrect";
}

// For the BPF, the coefficients are normalized so that the filter’s response at the mid‐frequency (f_mid) is 1.
// Compute the frequency response at f_mid and test that it is close to unity.
TEST_F(Filter_Fixture, BPF_Normalization) {
    int centre = (bpf_taps - 1) / 2;
    dy4::real f_mid = (f_low + f_high) / 2;
    dy4::real response = 0.0;
    // Compute the real part of the frequency response at f_mid (assuming linear phase)
    for (int k = 0; k < bpf_taps; k++) {
        int n = k - centre;
        response += bpf[k] * std::cos(2 * M_PI * n * f_mid / fs);
    }
    EXPECT_NEAR(response, 1.0, 1e-6)
        << "BPF normalization at f_mid is incorrect";
}

}  // namespace

// main() is provided by Google Test's framework via RUN_ALL_TESTS()
