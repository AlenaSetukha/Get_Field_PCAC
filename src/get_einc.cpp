#include <iostream>
#include <cmath>
#include <complex>

#include "element_geom.h"
#include "constants.h"

using namespace Constants;

//===============================================================================================
//---------------------Функция расчета падающего поля Einc(x) = e0 * e^i(k, x)-------------------
//===============================================================================================
void get_einc(const double (&x)[3], const double (&e0)[3], const double (&k_vec)[3], std::complex<double> (&e_inc)[3])
{
    std::complex<double> deg = pow(Constants::e, Constants::i_complex * scal_prod(k_vec, x));
    e_inc[0] = e0[0] * deg;
    e_inc[1] = e0[1] * deg;
    e_inc[2] = e0[2] * deg;
}