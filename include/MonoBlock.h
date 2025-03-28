#ifndef DY4_MONOBLOCK_H
#define DY4_MONOBLOCK_H

#include "fourier.h"
#include "filter.h"
#include "iofunc.h"
#include "dy4.h"

using namespace std;
using namespace dy4;

void upsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &upsampled);
void downsample(const std::vector<dy4::real> &signal, int factor, std::vector<dy4::real> &downsampled);
void generate_fm_values(std::vector<dy4::real>& signal, dy4::real min_val, dy4::real max_val);

// FIR filtering function
// void readIQData(const string& filename, std::vector<dy4::real>& iq_data);
// void firwin(int num_taps, dy4::real cutoff, dy4::real fs, std::vector<dy4::real>& coeffs);
// void processIQData(const std::vector<dy4::real>& iq_data, std::vector<dy4::real>& audio_data);
// void writeWAV(const string& filename, const std::vector<dy4::real>& audio_data, dy4::real sample_rate);
// void normalizeVector(std::std::vector<dy4::real> &vec, dy4::real new_min, dy4::real new_max);

#endif // DY4_MONOBLOCK_H