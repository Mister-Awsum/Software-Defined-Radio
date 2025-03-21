/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write psd unit tests, using Google C++ test framework.
// (it is based on https://github.com/google/googletest/blob/main/docs/index.md)

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "fourier.h"
#include "gtest/gtest.h"

namespace {

	class psd_Fixture: public ::testing::Test {

		public:

			const int N = 1024;
            const int Fs = 1e3;
			const int lower_bound = -1;
			const int upper_bound = 1;
			const real EPSILON = 1e-4;

			std::vector<real> x;
            
            std::vector<real> Freq_ref;
            std::vector<real> psdEst_ref;

            std::vector<real> Freq_test;
            std::vector<real> psdEst_test;

			// Twiddles will be used by the functions you will develop in-lab
			std::vector<std::complex<real>> Twiddle1D;
			std::vector<std::vector<std::complex<real>>> Twiddle2D;

			psd_Fixture( ) {
				x.resize(N);
				Twiddle1D.resize(N);
				Twiddle2D.resize(N, std::vector<std::complex<real>>(N));
			}

			void SetUp( ) {
				generate_random_values(x, lower_bound, upper_bound);
				estimatePSD(x, Fs, Freq_ref, psdEst_ref);
			}

			void TearDown( ) {
			}

			~psd_Fixture( )  {
			}
	};

    TEST_F(psd_Fixture, psd_matrix_NEAR) { 

		psd_matrix(x, Fs, Freq_test, psdEst_test);

		ASSERT_EQ(psdEst_ref.size(), psdEst_test.size()) << "The output vector sizes for psd_reference and psd_matrix are of unequal length";

		for (int i = 0; i < (int)psdEst_ref.size(); ++i) {
			EXPECT_NEAR(std::abs(psdEst_ref[i]), std::abs(psdEst_test[i]), EPSILON) << "Original/psd_matrix vectors differ at index " << i;
		}
	}

	TEST_F(psd_Fixture, psd_twiddle_2d_NEAR) { 

		psd_2d_twiddle(x, Fs, Freq_test, psdEst_test);

		ASSERT_EQ(psdEst_ref.size(), psdEst_test.size()) << "The output vector sizes for psd_reference and psd_2d_twiddle are of unequal length";

		for (int i = 0; i < (int)psdEst_ref.size(); ++i) {
			EXPECT_NEAR(std::abs(psdEst_ref[i]), std::abs(psdEst_test[i]), EPSILON) << "Original/psd_2d_twiddle vectors differ at index " << i;
		}
	}

	TEST_F(psd_Fixture, psd_unrolling_NEAR) { 

		psd_unrolling(x, Fs, Freq_test, psdEst_test);

		ASSERT_EQ(psdEst_ref.size(), psdEst_test.size()) << "The output vector sizes for psd_reference and psd_unrolling are of unequal length";

		for (int i = 0; i < (int)psdEst_ref.size(); ++i) {
			EXPECT_NEAR(std::abs(psdEst_ref[i]), std::abs(psdEst_test[i]), EPSILON) << "Original/psd_unrolling vectors differ at index " << i;
		}
	}


} // end of namespace

