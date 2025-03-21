/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write convolution unit tests, using Google C++ test framework.
// (it is based on https://github.com/google/googletest/blob/main/docs/index.md)

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"
#include "gtest/gtest.h"

namespace {

	class Convolution_Fixture: public ::testing::Test {

		public:

			const int N = 1024;	// signal size
			const int M = 101;	// kernel size
			const int lower_bound = -1;
			const int upper_bound = 1;
			const real EPSILON = 1e-4;

			std::vector<real> x, h, y_reference, y_test;

			Convolution_Fixture() {
				x.resize(N);
				h.resize(M);
				y_reference.resize(N + M - 1);
				y_test.resize(N + M - 1);
			}

			void SetUp() {
				generate_random_values(x, lower_bound, upper_bound);
				generate_random_values(h, lower_bound, upper_bound);
				convolveFIR_reference(y_reference, x, h);
			}

			void TearDown() {
			}

			~Convolution_Fixture() {
			}
	};

	TEST_F(Convolution_Fixture, convolveFIR_filterinefficient_NEAR) {

		convolveFIR_inefficient(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_inefficient are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_inefficient vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, convolveFIR_loop_x_NEAR) {

		convolveFIR_loop_x(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_loop_x are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_loop_x vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, convolveFIR_loop_h_NEAR) {

		convolveFIR_loop_h(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_loop_h are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_loop_h vectors differ at index " << i;
		}
	}
	
	TEST_F(Convolution_Fixture, convolveFIR_Valid_Range_NEAR) {

		convolveFIR_Valid_Range(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_Valid_Range are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_Valid_Range vectors differ at index " << i;
		}
	}
	
	TEST_F(Convolution_Fixture, convolveFIR_Range_and_Unrolling_NEAR) {

		convolveFIR_Range_and_Unrolling(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_Range_and_Unrolling are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_Range_and_Unrolling vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, convolveFIR_il4a_NEAR) {

		convolveFIR_il4a(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_il4a are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_il4a vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, convolveFIR_il4b_NEAR) {

		convolveFIR_il4b(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_il4b are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_il4b vectors differ at index " << i;
		}
	}

	TEST_F(Convolution_Fixture, convolveFIR_il4c_NEAR) {

		convolveFIR_il4c(y_test, x, h);

		ASSERT_EQ(y_reference.size(), y_test.size()) << "Output vector sizes for convolveFIR_reference and convolveFIR_il4c are unequal";

		for (int i = 0; i < (int)y_reference.size(); ++i) {
			EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) << "Original/convolveFIR_il4c vectors differ at index " << i;
		}
	}

} // end of namespace
