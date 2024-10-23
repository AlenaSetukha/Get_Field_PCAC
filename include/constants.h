#ifndef EM_FIELD_CONSTANTS_H
#define EM_FIELD_CONSTANTS_H

#include <complex>
#include <limits>
#include <iostream>

//========================================================================
//----------------------------Основные константы--------------------------
//========================================================================

namespace Constants {
    constexpr double ra = 57.295779513;
    constexpr double eps0 = 8.8541878128 * 1e-12;
    constexpr double m0 = 4 * M_PI * 1e-7;
    const double machine_zero = std::numeric_limits<double>::epsilon();
    constexpr int c_vacuum = 299792458; // м/с
} // namespace Constants
#endif // EM_FIELD_CONSTANTS_H
