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
    const dy4::real fs = 2.4e6;   // Example sample rate for audio (or IF)
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

// Test that the BPF coefficients match expected values
TEST_F(Filter_Fixture, BPF_Coefficients) {
    // Precomputed reference coefficients for the BPF
    std::vector<dy4::real> expected_bpf = {
        0.00000000e+00, 6.67282191e-07, -4.02388095e-05, -8.32372177e-05,
        -1.83865863e-05, -6.15626429e-19, -2.02139794e-04, -1.87651451e-04,
        5.24114541e-04, 1.20046597e-03, 6.58172315e-04, -4.99876493e-04,
        -4.69497700e-04, 3.55130316e-04, -5.82373837e-04, -3.13641941e-03,
        -3.21213063e-03, 7.43815632e-04, 3.87825834e-03, 2.21349595e-03,
        0.00000000e+00, 2.81339157e-03, 6.27244173e-03, 1.53434093e-03,
        -8.48116231e-03, -1.06521836e-02, -2.56059988e-03, 2.03801917e-03,
        -3.55275313e-03, -5.05132128e-03, 9.02236032e-03, 2.27664352e-02,
        1.40964271e-02, -7.38982159e-03, -1.21574976e-02, 4.13206737e-18,
        -3.14868754e-03, -2.86199747e-02, -3.52918219e-02, 2.66701372e-03,
        4.31924811e-02, 3.41664176e-02, 1.62267504e-03, 1.39972673e-02,
        5.72380741e-02, 2.78439567e-02, -1.04261330e-01, -1.97594944e-01,
        -9.58462984e-02, 1.40905388e-01, 2.66699896e-01, 1.40905388e-01,
        -9.58462984e-02, -1.97594944e-01, -1.04261330e-01, 2.78439567e-02,
        5.72380741e-02, 1.39972673e-02, 1.62267504e-03, 3.41664176e-02,
        4.31924811e-02, 2.66701372e-03, -3.52918219e-02, -2.86199747e-02,
        -3.14868754e-03, 4.13206737e-18, -1.21574976e-02, -7.38982159e-03,
        1.40964271e-02, 2.27664352e-02, 9.02236032e-03, -5.05132128e-03,
        -3.55275313e-03, 2.03801917e-03, -2.56059988e-03, -1.06521836e-02,
        -8.48116231e-03, 1.53434093e-03, 6.27244173e-03, 2.81339157e-03,
        0.00000000e+00, 2.21349595e-03, 3.87825834e-03, 7.43815632e-04,
        -3.21213063e-03, -3.13641941e-03, -5.82373837e-04, 3.55130316e-04,
        -4.69497700e-04, -4.99876493e-04, 6.58172315e-04, 1.20046597e-03,
        5.24114541e-04, -1.87651451e-04, -2.02139794e-04, -6.15626429e-19,
        -1.83865863e-05, -8.32372177e-05, -4.02388095e-05, 6.67282191e-07,
        0.00000000e+00
    };

    // Generate the BPF coefficients using the function under test
    std::vector<dy4::real> generated_bpf;
    generateBPF(generated_bpf, f_low, f_high, fs, bpf_taps);

    // Check that the generated coefficients match the expected coefficients
    ASSERT_EQ(generated_bpf.size(), expected_bpf.size())
        << "BPF coefficient size mismatch";

    for (size_t i = 0; i < expected_bpf.size(); i++) {
        EXPECT_NEAR(generated_bpf[i], expected_bpf[i], 1e-6)
            << "BPF coefficient mismatch at index " << i;
    }
}

}  // namespace

// main() is provided by Google Test's framework via RUN_ALL_TESTS()
