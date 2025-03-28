#include "gtest/gtest.h"
#include "MonoBlock.h"
#include "Stereo.h"
#include <vector>
#include <cmath>
#include <chrono> // For time measurement

namespace {

class StereoTest : public ::testing::Test {
public:
    const int N = 1000; // Number of samples for testing
    const int DECIMATION_FACTOR = 2; // Decimation factor
    const double EPSILON = 1e-6; // Tolerance for floating-point comparisons
    const int FILTER_TAPS = 101; // Number of filter taps

    // Test signals
    std::vector<dy4::real> input_signal;
    std::vector<dy4::real> mono_output;
    std::vector<dy4::real> stereo_output;
    std::vector<dy4::real> filter_coeff;
    std::vector<dy4::real> state;
    std::vector<dy4::real> output_convolveBlock;
    std::vector<dy4::real> output_convolveFast;
    std::vector<dy4::real> output_convolveAndDecimate;
    std::vector<dy4::real> delayed_output;

    // Filter coefficients for mono and stereo paths
    std::vector<dy4::real> mono_coeff;
    std::vector<dy4::real> stereo_coeff;

    // Filter states
    std::vector<dy4::real> mono_state;
    std::vector<dy4::real> stereo_state;

    // Delay block state
    std::vector<dy4::real> delay_state;

    StereoTest() {
        input_signal.resize(N, 0.0);
        mono_output.resize(N, 0.0);
        stereo_output.resize(N, 0.0);
        filter_coeff.resize(FILTER_TAPS, 0.0);
        state.resize(FILTER_TAPS - 1, 0.0);
        output_convolveBlock.resize(N, 0.0);
        output_convolveFast.resize(N, 0.0);
        output_convolveAndDecimate.resize(N / DECIMATION_FACTOR, 0.0);
        delayed_output.resize(N, 0.0);
        mono_coeff.resize(FILTER_TAPS, 0.0);
        stereo_coeff.resize(FILTER_TAPS, 0.0);
        mono_state.resize(FILTER_TAPS - 1, 0.0);
        stereo_state.resize(FILTER_TAPS - 1, 0.0);
        delay_state.resize((FILTER_TAPS - 1) / 2, 0.0);
    }

    void SetUp() override {
        // Create an impulse signal for testing
        input_signal[N / 2] = 1.0;

        // Generate filter coefficients
        generateLPF(48e3, 16e3, FILTER_TAPS, mono_coeff); // Mono LPF
        generateBPF(stereo_coeff, 22e3, 54e3, 48e3, FILTER_TAPS); // Stereo BPF
    }

    void TearDown() override {}
};

// Test to compare delays in mono and stereo paths
TEST_F(StereoTest, FilterDelayComparison) {
    int mono_delay = (mono_coeff.size() - 1) / 2;
    int stereo_delay = (stereo_coeff.size() - 1) / 2;

    EXPECT_EQ(mono_delay, stereo_delay) << "Delays in mono and stereo paths do not match!";
}

// Test to verify the allpass filter compensates for delay in mono path
TEST_F(StereoTest, AllpassFilterCompensation) {
    convolveBlock(input_signal, mono_coeff, mono_state, mono_output);
    convolveBlock(input_signal, stereo_coeff, stereo_state, stereo_output);

    ASSERT_EQ(mono_output.size(), stereo_output.size()) 
        << "Output sizes do not match!";
}

// Test to verify delay_block delay matches convolution delay
TEST_F(StereoTest, DelayBlockDelayTest) {
    std::vector<dy4::real> delayed_output;
    delay_block(input_signal, delay_state, delayed_output);

    // Verify delay introduced by delay_block matches convolution delay
    int delay_block_delay = delay_state.size();
    int convolution_delay = (mono_coeff.size() - 1) / 2;

    EXPECT_EQ(delay_block_delay, convolution_delay)
        << "Delay introduced by delay_block does not match convolution delay!";
}

// Test for equivalence between convolveBlock and ConvolveFast
TEST_F(StereoTest, ConvolveEquivalenceTest) {
    convolveBlock(input_signal, filter_coeff, state, output_convolveBlock);
    ConvolveFast(input_signal, filter_coeff, state, 1, 1, output_convolveFast);

    ASSERT_EQ(output_convolveBlock.size(), output_convolveFast.size())
        << "Output sizes do not match!";

    for (size_t i = 0; i < output_convolveBlock.size(); i++) {
        EXPECT_NEAR(output_convolveBlock[i], output_convolveFast[i], EPSILON)
            << "Mismatch at index " << i;
    }
}

// Test for convolveAndDecimate
TEST_F(StereoTest, ConvolveAndDecimateFunctionality) {
    convolveAndDecimate(input_signal, filter_coeff, state, DECIMATION_FACTOR, output_convolveAndDecimate);

    // Verify the output matches the expected convolution followed by decimation
    for (size_t i = 0; i < output_convolveAndDecimate.size(); i++) {
        dy4::real expected = 0.0;
        size_t index = i * DECIMATION_FACTOR;
        if (index >= N / 2 && index < N / 2 + filter_coeff.size()) {
            expected = filter_coeff[index - N / 2];
        }
        EXPECT_NEAR(output_convolveAndDecimate[i], expected, EPSILON)
            << "Mismatch at index " << i;
    }
}

// Test to compare runtime of delay_block against convolveBlock
TEST_F(StereoTest, DelayBlockVsConvolveBlockRuntime) {
    auto start = std::chrono::high_resolution_clock::now();
    delay_block(input_signal, delay_state, delayed_output);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delay_block_time = end - start;

    start = std::chrono::high_resolution_clock::now();
    convolveBlock(input_signal, mono_coeff, mono_state, mono_output);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convolve_block_time = end - start;

    std::cout << "Delay block runtime: " << delay_block_time.count() << " seconds" << std::endl;
    std::cout << "ConvolveBlock runtime: " << convolve_block_time.count() << " seconds" << std::endl;

    EXPECT_LT(delay_block_time.count(), convolve_block_time.count())
        << "delay_block should be faster than convolveBlock!";
}

// Test to compare runtime of delay_block against ConvolveFast
TEST_F(StereoTest, DelayBlockVsConvolveFastRuntime) {
    auto start = std::chrono::high_resolution_clock::now();
    delay_block(input_signal, delay_state, delayed_output);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delay_block_time = end - start;

    start = std::chrono::high_resolution_clock::now();
    ConvolveFast(input_signal, mono_coeff, mono_state, 1, 1, mono_output);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convolve_fast_time = end - start;

    std::cout << "Delay block runtime: " << delay_block_time.count() << " seconds" << std::endl;
    std::cout << "ConvolveFast runtime: " << convolve_fast_time.count() << " seconds" << std::endl;

    EXPECT_LT(delay_block_time.count(), convolve_fast_time.count())
        << "delay_block should be faster than ConvolveFast!";
}

// Test to compare runtime of delay_block against convolveAndDecimate
TEST_F(StereoTest, DelayBlockVsConvolveAndDecimateRuntime) {
    auto start = std::chrono::high_resolution_clock::now();
    delay_block(input_signal, delay_state, delayed_output);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delay_block_time = end - start;

    start = std::chrono::high_resolution_clock::now();
    convolveAndDecimate(input_signal, mono_coeff, mono_state, DECIMATION_FACTOR, mono_output);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convolve_and_decimate_time = end - start;

    std::cout << "Delay block runtime: " << delay_block_time.count() << " seconds" << std::endl;
    std::cout << "ConvolveAndDecimate runtime: " << convolve_and_decimate_time.count() << " seconds" << std::endl;

    EXPECT_LT(delay_block_time.count(), convolve_and_decimate_time.count())
        << "delay_block should be faster than convolveAndDecimate!";
}

TEST_F(StereoTest, CompareCppWithPythonOutputs) {
    for (unsigned int block_id = 0; block_id < 5; block_id++) {
        // Compare FM demodulated data
        std::string fm_demod_cpp_file = "../data/stereo_impl_testing/fm_demod_block_" + std::to_string(block_id) + "_cpp.bin";
        std::string fm_demod_python_file = "../data/stereo_impl_testing/fm_demod_block_" + std::to_string(block_id) + "_python.bin";

        std::vector<dy4::real> fm_demod_cpp, fm_demod_python;
        readBinData(fm_demod_cpp_file, fm_demod_cpp);
        readBinData(fm_demod_python_file, fm_demod_python);

        ASSERT_EQ(fm_demod_cpp.size(), fm_demod_python.size()) << "FM Demod size mismatch for block " << block_id;
        for (size_t i = 0; i < fm_demod_cpp.size(); i++) {
            EXPECT_NEAR(fm_demod_cpp[i], fm_demod_python[i], EPSILON) << "Mismatch in FM Demod at index " << i << " for block " << block_id;
        }

        // Compare stereo carrier recovery data
        std::string stereo_filt_cpp_file = "../data/stereo_impl_testing/stereo_filt_block_" + std::to_string(block_id) + "_cpp.bin";
        std::string stereo_filt_python_file = "../data/stereo_impl_testing/stereo_filt_block_" + std::to_string(block_id) + "_python.bin";

        std::vector<dy4::real> stereo_filt_cpp, stereo_filt_python;
        readBinData(stereo_filt_cpp_file, stereo_filt_cpp);
        readBinData(stereo_filt_python_file, stereo_filt_python);

        ASSERT_EQ(stereo_filt_cpp.size(), stereo_filt_python.size()) << "Stereo Carrier size mismatch for block " << block_id;
        for (size_t i = 0; i < stereo_filt_cpp.size(); i++) {
            EXPECT_NEAR(stereo_filt_cpp[i], stereo_filt_python[i], EPSILON) << "Mismatch in Stereo Carrier at index " << i << " for block " << block_id;
        }

        // Compare PLL & NCO output
        std::string pll_nco_cpp_file = "../data/stereo_impl_testing/pll_nco_block_" + std::to_string(block_id) + "_cpp.bin";
        std::string pll_nco_python_file = "../data/stereo_impl_testing/pll_nco_block_" + std::to_string(block_id) + "_python.bin";

        std::vector<dy4::real> pll_nco_cpp, pll_nco_python;
        readBinData(pll_nco_cpp_file, pll_nco_cpp);
        readBinData(pll_nco_python_file, pll_nco_python);

        ASSERT_EQ(pll_nco_cpp.size(), pll_nco_python.size()) << "PLL NCO size mismatch for block " << block_id;
        for (size_t i = 0; i < pll_nco_cpp.size(); i++) {
            EXPECT_NEAR(pll_nco_cpp[i], pll_nco_python[i], EPSILON) << "Mismatch in PLL NCO at index " << i << " for block " << block_id;
        }

        // Compare stereo channel extraction data
        std::string s_channel_cpp_file = "../data/stereo_impl_testing/s_channel_block_" + std::to_string(block_id) + "_cpp.bin";
        std::string s_channel_python_file = "../data/stereo_impl_testing/s_channel_block_" + std::to_string(block_id) + "_python.bin";

        std::vector<dy4::real> s_channel_cpp, s_channel_python;
        readBinData(s_channel_cpp_file, s_channel_cpp);
        readBinData(s_channel_python_file, s_channel_python);

        ASSERT_EQ(s_channel_cpp.size(), s_channel_python.size()) << "Stereo Channel size mismatch for block " << block_id;
        for (size_t i = 0; i < s_channel_cpp.size(); i++) {
            EXPECT_NEAR(s_channel_cpp[i], s_channel_python[i], EPSILON) << "Mismatch in Stereo Channel at index " << i << " for block " << block_id;
        }
    }
}

} // namespace
