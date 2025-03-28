#ifndef DY4_RFFRONTEND_H
#define DY4_RFFRONTEND_H

#include "fourier.h"
#include "filter.h"
#include "iofunc.h"
#include "dy4.h"
#include "SyncQueue.h"

using namespace std;
using namespace dy4;

void fmDemodulate(const std::vector<dy4::real>& i_ds, const std::vector<dy4::real>& q_ds, std::vector<dy4::real>& state_phase, std::vector<dy4::real>& fm_demod);
void fmDemodArctan(const std::vector<dy4::real>& i_ds, const std::vector<dy4::real>& q_ds, dy4::real& state_phase, std::vector<dy4::real>& fm_demod);
void RFProcessingThread(SyncQueue &syncQueue, const std::vector<dy4::real> &rf_coeff, int if_decim,
                        unsigned int u_read_size, int rf_taps);
void splitIQSamples(const std::vector<dy4::real> &block_data, std::vector<dy4::real> &i_samples, std::vector<dy4::real> &, unsigned int u_read_size);
// Read data from stdin
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<dy4::real> &block_data);

#endif // DY4_RFRONTEND_H
