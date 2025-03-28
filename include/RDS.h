#ifndef RDS_H
#define RDS_H

#include <vector>
#include "dy4.h"
#include "filter.h"

void extract_rds_channel(const std::vector<dy4::real> &fm_demod,
                         const std::vector<dy4::real> &rds_coeff,
                         std::vector<dy4::real> &state_rds,
                         std::vector<dy4::real> &rds_filtered);

void extract_rds_carrier(const std::vector<dy4::real> &rds_filtered,
                         const std::vector<dy4::real> &carrier_coeff,
                         std::vector<dy4::real> &state_carrier,
                         std::vector<dy4::real> &carrier);

void downconvert_rds_to_baseband(const std::vector<dy4::real> &rds_filtered,
                                  const std::vector<dy4::real> &carrier,
                                  std::vector<dy4::real> &baseband);

#endif