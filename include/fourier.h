/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FOURIER_H
#define DY4_FOURIER_H

// Add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "dy4.h"

using namespace dy4;

// Declaration of function prototypes
void DFT(const std::vector<real> &,
	std::vector<std::complex<real>> &);

// You should add your own IDFT
// time-permitting you can build your own function for FFT

void computeVectorMagnitude(const std::vector<std::complex<real>> &,
	std::vector<real> &);

// Provide the prototype to estimate PSD
void estimatePSD(const std::vector<real> &samples, unsigned int Fs, std::vector<real> &freq, std::vector<real> &psd_est);

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

// Prototype DFT optimizations functions
void DFT_reference(const std::vector<real>& x, std::vector<std::complex<real>>& Xf);
void DFT_init_bins(const std::vector<real>& x, std::vector<std::complex<real>>& Xf);
void DFT_code_motion(const std::vector<real> &x, std::vector<std::complex<real>> &Xf);
void DFT_1d_twiddle(const std::vector<real> &x, std::vector<std::complex<real>> &Twiddle1D, std::vector<std::complex<real>> &Xf);
void DFT_2d_twiddle(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf);
void DFT_loop_m(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf);
void DFT_loop_k(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf);
void DFT_unrolling(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf);

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D);
void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>>& Twiddle2D);

// Prototype matrix PSD functions
void psd_matrix(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst);
void psd_2d_twiddle(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst);
void psd_unrolling(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst);

#endif // DY4_FOURIER_H
