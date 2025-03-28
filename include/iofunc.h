/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// Add headers as needed
#include <iomanip>
#include "dy4.h"

// Declaration of function prototypes
void printRealVector(const std::vector<dy4::real> &);

void printComplexVector(const std::vector<std::complex<dy4::real>> &);

void readBinData(const std::string, std::vector<dy4::real> &);

void writeBinData(const std::string, const std::vector<dy4::real> &);

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void generate_random_values(std::vector<dy4::real>& x, const dy4::real& lower_bound, const dy4::real& upper_bound);
void appendBinData(const std::string out_fname, const std::vector<dy4::real> &bin_data);

#endif // DY4_IOFUNC_H
