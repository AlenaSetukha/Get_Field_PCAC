#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

namespace Constants
{
    constexpr double pi = 3.1415926535;                       // pi
    constexpr double pi_reverse =  0.07957747154;             // 1/4pi
    constexpr double ra = 57.295779513;
    constexpr double eps0 = 8.8541878128 * 1e-12;
    constexpr double m0 = 4 * pi * 1e-7;
    constexpr std::complex<double> i_complex = std::complex<double>(0.0, 1.0);
    constexpr double e = 2.7182818284;
    constexpr double machine_zero = 1e-16;
}
#endif
