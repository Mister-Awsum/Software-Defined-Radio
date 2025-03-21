/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_GENFUNC_H
#define DY4_GENFUNC_H

// Add headers as needed
#include <iostream>
#include <vector>
#include <complex>

#include "dy4.h"

using namespace dy4;

// Declaration of function prototypes
void generateSin(std::vector<dy4::real> &, std::vector<dy4::real> &, dy4::real, \
	real, real, real, real);

void addSin(const std::vector<std::vector<dy4::real>> &, std::vector<dy4::real> &);

void generateRandomSamples(std::vector<dy4::real> &, unsigned int, \
	unsigned short int, unsigned char);



#endif // DY4_GENFUNC_H
