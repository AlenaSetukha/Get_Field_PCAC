#include <iostream>
#include <cmath>
#include <complex>

#include "element_geom.h"
#include "constants.h"

//===============================================================================================
//---------------------Функция расчета падающего поля Einc(x) = e0 * e^i(k, x)-------------------
//===============================================================================================

void get_einc(const double* x, const double* e0, const double* k_vec, std::complex<double>* e_inc)
{
    Constants c;
    std::complex<double> deg = pow(c.e, c.i_complex * scal_prod(k_vec, x));
    e_inc[0] = e0[0] * deg;
    e_inc[1] = e0[1] * deg;
    e_inc[2] = e0[2] * deg;
}