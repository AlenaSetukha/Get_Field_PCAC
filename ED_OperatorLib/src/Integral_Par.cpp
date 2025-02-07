#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "Integral_Par.h"

Integral_Par::Integral_Par(const int idim_in, const int n_start_in, const int p_max_in, const double eps_in) {
    idim = idim_in;
    n_start = n_start_in;
    p_max = p_max_in;
    eps = eps_in;
}

Integral_Par::Integral_Par(const Integral_Par& obj) {
    idim = obj.idim;
    n_start = obj.n_start;
    p_max = obj.p_max;
    eps = obj.eps;
}

