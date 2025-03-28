/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write DFT benchmark functions, using Google C++ benchmark framework.
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

static void Bench_DFT_reference(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		DFT_reference(x, Xf);
	}
}

// register benchmark Bench_DFT_reference //

BENCHMARK(Bench_DFT_reference)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_init_bins(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		DFT_init_bins(x, Xf);
	}
}

// register benchmark Bench_DFT_init_bins //

BENCHMARK(Bench_DFT_init_bins)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_code_motion(benchmark::State& state) {
	int N = state.range(0);

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);

	for (auto _ : state) {
		DFT_code_motion(x, Xf);
	}
}

// register benchmark Bench_DFT_code_motion //

BENCHMARK(Bench_DFT_code_motion)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_1d_twiddle(benchmark::State& state) {
	int N = state.range(0);

	std::vector<std::complex<real>> Twiddle1D;
	Twiddle1D.resize(N);

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);
	generate_DFT_twiddles(N, Twiddle1D);

	for (auto _ : state) {
		DFT_1d_twiddle(x, Twiddle1D, Xf);
	}
}

// register benchmark Bench_DFT_1d_twiddle //

BENCHMARK(Bench_DFT_1d_twiddle)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_2d_twiddle(benchmark::State& state) {
	int N = state.range(0);

	std::vector<std::vector<std::complex<real>>> Twiddle2D;
	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);
	generate_DFT_matrix(N, Twiddle2D);

	for (auto _ : state) {
		DFT_2d_twiddle(x, Twiddle2D, Xf);
	}
}

// register benchmark Bench_DFT_2d_twiddle //

BENCHMARK(Bench_DFT_2d_twiddle)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_loop_m(benchmark::State& state) {
	int N = state.range(0);

	std::vector<std::vector<std::complex<real>>> Twiddle2D;
	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);
	generate_DFT_matrix(N, Twiddle2D);

	for (auto _ : state) {
		DFT_loop_m(x, Twiddle2D, Xf);
	}
}

// register benchmark Bench_DFT_loop_m //

BENCHMARK(Bench_DFT_loop_m)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_loop_k(benchmark::State& state) {
	int N = state.range(0);

	std::vector<std::vector<std::complex<real>>> Twiddle2D;
	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);
	generate_DFT_matrix(N, Twiddle2D);

	for (auto _ : state) {
		DFT_loop_k(x, Twiddle2D, Xf);
	}
}

// register benchmark Bench_DFT_loop_k //

BENCHMARK(Bench_DFT_loop_k)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

static void Bench_DFT_unrolling(benchmark::State& state) {
	int N = state.range(0);

	std::vector<std::vector<std::complex<real>>> Twiddle2D;
	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));

	std::vector<real> x(N);
	std::vector<std::complex<real>> Xf(N);
	generate_random_values(x, lower_bound, upper_bound);
	generate_DFT_matrix(N, Twiddle2D);

	for (auto _ : state) {
		DFT_unrolling(x, Twiddle2D, Xf);
	}
}

// register benchmark Bench_DFT_unrolling //

BENCHMARK(Bench_DFT_unrolling)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

////////////////////////////////////////////

