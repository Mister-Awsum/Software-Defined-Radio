#ifndef DY4_STEREO_H
#define DY4_STEREO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include <cstdint>

#include "fourier.h"
#include "filter.h"
#include "iofunc.h"
#include "MonoBlock.h"
#include "dy4.h"

using namespace std;
using namespace dy4;

const dy4::real normBandwidth = 0.01;   // Normalized bandwidth for the PLL

// PLL and NCO functionality
void fmPLL(const std::vector<dy4::real>& pllIn, dy4::real freq, dy4::real Fs, dy4::real ncoScale, dy4::real phaseAdjust, dy4::real normBandwidth, std::vector<dy4::real>& output);

// Stereo Carrier Recovery Function
void s_carrier_recovery(vector<dy4::real> i_sig, vector<dy4::real> coeff, vector<dy4::real> state, dy4::real Fs, vector<dy4::real>& output);
void s_processing(vector<dy4::real> s_carrier, vector<dy4::real> s_channel, vector<short int> &audio_data);

#endif // DY4_STEREO_H