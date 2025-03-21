/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_LOGFUNC_H
#define DY4_LOGFUNC_H

// Add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "dy4.h"

void genIndexVector(std::vector<dy4::real> &, \
	const int);

void logVector(const std::string, \
	const std::vector<dy4::real> &, \
	const std::vector<dy4::real> &);

#endif // DY4_LOGFUNC_H
