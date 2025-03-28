#ifndef DY4_AUDIOPROCESSING_H
#define DY4_AUDIOPROCESSING_H

#include "SyncQueue.h"
#include "Stereo.h"
#include "MonoBlock.h"
#include "filter.h"
#include "dy4.h"

void AudioProcessingThread(SyncQueue &syncQueue, const std::vector<dy4::real> &if_coeff, 
                           std::vector<dy4::real> &state_audio, int U, int D, 
                           const std::vector<dy4::real> &extract_coeff, std::vector<dy4::real> &state_extract, std::vector<dy4::real> &stereo_process_state,
                           const std::vector<dy4::real> &pilot_coeff, std::vector<dy4::real> &state_recovery, std::vector<dy4::real> &state_block, std::vector<dy4::real> &output_block,
                           dy4::real pll_freq, dy4::real pll_phase, dy4::real if_Fs, dy4::real normBandwidth, dy4::real nco_Scalar,
                           std::vector<dy4::real> &pll_states, const std::string &mode);

#endif // DY4_AUDIOPROCESSING_H