#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

class Constants
{
public:
    const double pi = 3.1415926535;                       // pi
    const double pi_reverse =  0.07957747154;             // 1/4pi
    const double ra = 57.295779513;
    const double eps0 = 8.8541878128 * pow(10.0, -12.0);
    const double m0 = 4 * pi * pow(10.0, -7.0);
    const std::complex<double> i_complex = std::complex<double>(0.0, 1.0);
    const double e = 2.7182818284;
    const double machine_zero = pow(10.0, -16);
};
#endif
