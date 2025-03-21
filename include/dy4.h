/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_DY4_H
#define DY4_DY4_H

#include <cmath>  // For math functions
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
#include <cstdint>
#include <unistd.h>

// Macro for conditional compilation
namespace dy4 {
    #ifdef DOUBLE
        typedef double real;
    #else
        typedef float real;
    #endif
}

// General constants
const dy4::real PI = 3.14159265358979323846;
const int NFFT = 512;  // Number of FFT points

#endif // DY4_DY4_H
