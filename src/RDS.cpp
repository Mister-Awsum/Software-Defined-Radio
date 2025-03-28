#include "RDS.h"
#include <cmath>
#include <iostream>

// Step 1: Filter FM-demod signal for 54–60 kHz band
void extract_rds_channel(const std::vector<dy4::real> &fm_demod,
                         const std::vector<dy4::real> &rds_coeff,
                         std::vector<dy4::real> &state_rds,
                         std::vector<dy4::real> &rds_filtered) {
    convolveAndDecimate(fm_demod, rds_coeff, state_rds, 1, rds_filtered);
}

// Step 2: Extract carrier ~57 kHz using bandpass filter
void extract_rds_carrier(const std::vector<dy4::real> &rds_filtered,
                         const std::vector<dy4::real> &carrier_coeff,
                         std::vector<dy4::real> &state_carrier,
                         std::vector<dy4::real> &carrier) {
    convolveAndDecimate(rds_filtered, carrier_coeff, state_carrier, 1, carrier);
}

// Step 3: Multiply to baseband
void downconvert_rds_to_baseband(const std::vector<dy4::real> &rds_filtered,
                                  const std::vector<dy4::real> &carrier,
                                  std::vector<dy4::real> &baseband) {
    baseband.resize(rds_filtered.size());
    for (size_t i = 0; i < rds_filtered.size(); i++) {
        baseband[i] = rds_filtered[i] * carrier[i];  // Mix down to baseband
    }
}