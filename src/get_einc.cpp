#include <cmath>
#include <complex>
#include <iostream>

#include "element_geom.h"

//===============================================================================================
//---------------------Функция расчета падающего поля Einc(x) = e0 * e^i(k, x)-------------------
//===============================================================================================
void get_einc(const double *x, const double *e0, const double *k_vec,
              std::complex<double> *e_inc) {
  std::complex<double> deg =
      std::exp(std::complex<double>(0., 1.) * scal_prod(k_vec, x));
  e_inc[0] = e0[0] * deg;
  e_inc[1] = e0[1] * deg;
  e_inc[2] = e0[2] * deg;
}
