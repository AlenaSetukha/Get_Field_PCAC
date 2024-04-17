#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "f_par.h"

f_par::f_par(const double rs_in)
{
    rs = rs_in;
}

f_par::f_par(const double rs_in, const std::complex<double> k_in)
{
    rs = rs_in;
    k = k_in;
}

f_par::f_par(const double rs_in, const std::complex<double> k_in, const double (&n_in)[3])
{
    rs = rs_in;
    k = k_in;
    n[0] = n_in[0];
    n[1] = n_in[1];
    n[2] = n_in[2];
}
