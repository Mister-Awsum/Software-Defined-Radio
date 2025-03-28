/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write psd benchmark functions, using Google C++ benchmark framework.
// (it is based on https://github.com/google/benchmark/blob/main/docs/user_guide.md)

#include <benchmark/benchmark.h>
#include "utils_bench.h"
#include "dy4.h"
#include "iofunc.h"
#include "fourier.h"

#define RANGE_MULTIPLIER 2
#define MIN_INPUT_SIZE 256
#define MAX_INPUT_SIZE (8 * MIN_INPUT_SIZE)

const int lower_bound = -1;
const int upper_bound = 1;
const int N = 1024;
const int Fs = 1e3;

// Reference Function

static void Bench_psd_reference(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
   std::vector<real> Freq_ref;
   std::vector<real> psdEst_ref;
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		estimatePSD(x, Fs, Freq_ref, psdEst_ref);
	}
}

// register benchmark Bench_psd_reference //

BENCHMARK(Bench_psd_reference)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_psd_matrix(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
   std::vector<real> Freq_test;
   std::vector<real> psdEst_test;
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		psd_matrix(x, Fs, Freq_test, psdEst_test);
	}
}

// register benchmark Bench_psd_matrix //

BENCHMARK(Bench_psd_matrix)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

//////////////////////////////////////////

static void Bench_psd_2d_twiddle(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
   std::vector<real> Freq_test;
   std::vector<real> psdEst_test;
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		psd_2d_twiddle(x, Fs, Freq_test, psdEst_test);
	}
}

// register benchmark Bench_psd_2d_twiddle //

BENCHMARK(Bench_psd_2d_twiddle)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_psd_unrolling(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
   std::vector<real> Freq_test;
   std::vector<real> psdEst_test;
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		psd_unrolling(x, Fs, Freq_test, psdEst_test);
	}
}

// register benchmark Bench_psd_unrolling //

BENCHMARK(Bench_psd_unrolling)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);