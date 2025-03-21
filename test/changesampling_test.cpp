/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include "dy4.h"           // defines dy4::real (e.g., float or double)
#include "gtest/gtest.h"
#include "MonoBlock.h"
#include <vector>
#include <cmath>
#include <algorithm>

// Assuming the upsample and downsample functions are declared as follows:
void upsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &upsampled);
void downsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &downsampled);

namespace {

class Resample_Fixture : public ::testing::Test {
public:
    // Use a short signal for testing
    const int N = 100;  
    const double EPSILON = 1e-6;
    std::vector<dy4::real> signal;

    Resample_Fixture() {
        signal.resize(N);
    }

    void SetUp() override {
        // Fill the signal with known values (here, a simple increasing sequence)
        for (int i = 0; i < N; i++) {
            signal[i] = static_cast<dy4::real>(i);  
        }
    }

    void TearDown() override {}
};

TEST_F(Resample_Fixture, UpsampleTest) {
    std::vector<dy4::real> upsampled_signal;
    int factor = 3;
    
    // Call updated upsample function
    upsample(signal, factor, upsampled_signal);

    // Check that the new size is N * factor
    EXPECT_EQ(upsampled_signal.size(), signal.size() * static_cast<size_t>(factor))
        << "Upsampled signal size is incorrect.";

    // Check that every factor-th element matches the original, and the inserted ones are zeros
    for (size_t i = 0; i < signal.size(); i++) {
        EXPECT_NEAR(upsampled_signal[i * factor], signal[i], EPSILON)
            << "Mismatch at index " << i * factor << ".";
        for (int j = 1; j < factor; j++) {
            EXPECT_NEAR(upsampled_signal[i * factor + j], 0.0, EPSILON)
                << "Expected zero at index " << (i * factor + j) << ".";
        }
    }
}

TEST_F(Resample_Fixture, DownsampleTest) {
    int factor = 4;
    std::vector<dy4::real> downsampled_signal;
    
    // Call updated downsample function
    downsample(signal, factor, downsampled_signal);

    // For an array of size N divisible by factor, the new size should be N / factor.
    size_t expected_size = signal.size() / factor;
    EXPECT_EQ(downsampled_signal.size(), expected_size)
        << "Downsampled signal size is incorrect.";

    // Check that the output contains every factor-th element of the original.
    for (size_t i = 0; i < expected_size; i++) {
        EXPECT_NEAR(downsampled_signal[i], signal[i * factor], EPSILON)
            << "Mismatch at downsampled index " << i << ".";
    }
}

TEST_F(Resample_Fixture, UpsampleDownsampleCombination) {
    int factor = 3;
    std::vector<dy4::real> upsampled_signal, downsampled_signal;
    
    // Upsample the signal
    upsample(signal, factor, upsampled_signal);
    // Then downsample by the same factor
    downsample(upsampled_signal, factor, downsampled_signal);

    // The size should return to original size.
    EXPECT_EQ(downsampled_signal.size(), signal.size())
        << "Combined upsample and downsample did not restore original size.";

    // Check that the values are (approximately) equal.
    for (size_t i = 0; i < signal.size(); i++) {
        EXPECT_NEAR(downsampled_signal[i], signal[i], EPSILON)
            << "Mismatch at index " << i << " after combined resampling.";
    }
}

}  // namespace

// main() is provided by Google Test's framework via RUN_ALL_TESTS().
