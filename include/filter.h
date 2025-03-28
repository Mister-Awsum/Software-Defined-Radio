/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// Add headers as needed
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "dy4.h"

using namespace dy4;

// Declaration of function prototypes
void generateLPF(real, real, unsigned short int, std::vector<real> &);
void generateHPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h);
void generateBPF(std::vector<dy4::real> &h, dy4::real f_low, dy4::real f_high, dy4::real fs, int taps);
void convolveFIR(std::vector<real> &, const std::vector<real> &, const std::vector<real> &);
void convolveBlock(const std::vector<dy4::real>& signal, const std::vector<dy4::real>& coeffs, std::vector<dy4::real>& state, std::vector<dy4::real>& output);
void convolveAndDecimate(const std::vector<dy4::real>& signal, const std::vector<dy4::real>& coeffs, std::vector<dy4::real>& state, int decimation_factor, std::vector<dy4::real>& output);
void ConvolveFast(const std::vector<dy4::real> &x, const std::vector<dy4::real> &h, std::vector<dy4::real> &state, const int &u, const int &d, std::vector<dy4::real> &y);

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);

void convolveFIR_loop_x(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_loop_h(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_Valid_Range(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_Range_and_Unrolling(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_il4a(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_il4b(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_il4c(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);


#endif // DY4_FILTER_H

