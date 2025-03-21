/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include <vector>
#include <cmath>
#include <iostream>
#include "gtest/gtest.h"
#include "dy4.h"         // Defines dy4::real (e.g., typedef double)
#include "MonoBlock.h"   // Contains declarations for convolveBlock, convolveAndDecimate, and ConvolveFast

// Helper function to compare two vectors with a given tolerance.
void ExpectVectorsNear(const std::vector<dy4::real>& expected,
                       const std::vector<dy4::real>& actual,
                       dy4::real tol) {
    ASSERT_EQ(expected.size(), actual.size()) << "Vector sizes differ";
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], actual[i], tol) << "Mismatch at index " << i;
    }
}

namespace {

class ConvolveFixture : public ::testing::Test {
protected:
    // Our test signal, filter coefficients, and state.
    // Signal: [1, 2, 3, 4, 5]
    // Coeffs: [0.2, 0.3, 0.5] (filter length 3, so state should have 2 samples)
    // Initial state: [0.1, 0.2]
    std::vector<dy4::real> signal;
    std::vector<dy4::real> coeffs;
    std::vector<dy4::real> state;
    std::vector<dy4::real> output;
    
    // Expected output computed manually:
    // [0.31, 0.8, 1.7, 2.7, 3.7]
    std::vector<dy4::real> expectedOutput;
    // Expected updated state should be the last (coeffs.size()-1) = 2 samples of signal: [4, 5]
    std::vector<dy4::real> expectedState;
    
    ConvolveFixture() {
        signal = {1, 2, 3, 4, 5};
        coeffs = {0.2, 0.3, 0.5};
        state = {0.1, 0.2};  // initial state
        
        expectedOutput = {0.31, 0.8, 1.7, 2.7, 3.7};
        expectedState = {4, 5};
    }
};

TEST_F(ConvolveFixture, ConvolveBlockProducesExpectedOutput) {
    // Preallocate output vector to the same size as signal.
    output.resize(signal.size(), 0.0);
    convolveBlock(signal, coeffs, state, output);
    ExpectVectorsNear(expectedOutput, output, 1e-6);
}

TEST_F(ConvolveFixture, ConvolveBlockUpdatesStateCorrectly) {
    // Reset state to initial value before test.
    state = {0.1, 0.2};
    output.resize(signal.size(), 0.0);
    convolveBlock(signal, coeffs, state, output);
    ExpectVectorsNear(expectedState, state, 1e-6);
}

TEST_F(ConvolveFixture, ConvolveAndDecimateProducesExpectedOutput) {
    // Use decimation factor d = 1 (i.e., no decimation) to compare with convolveBlock.
    int d = 1;
    output.resize(signal.size() / d, 0.0);
    // Reset state for a fair test.
    state = {0.1, 0.2};
    convolveAndDecimate(signal, coeffs, state, d, output);
    ExpectVectorsNear(expectedOutput, output, 1e-6);
    ExpectVectorsNear(expectedState, state, 1e-6);
}

TEST_F(ConvolveFixture, ConvolveFastProducesExpectedOutputFor_u1_d1) {
    // Using u = 1 and d = 1 should ideally give the same result as convolveBlock.
    int u = 1, d = 1;
    output.resize(signal.size(), 0.0);
    // Reset state
    state = {0.1, 0.2};
    ConvolveFast(signal, coeffs, state, u, d, output);
    ExpectVectorsNear(expectedOutput, output, 1e-6);
    ExpectVectorsNear(expectedState, state, 1e-6);
}

}  // namespace

// main() is provided by Google Test's framework via RUN_ALL_TESTS().
