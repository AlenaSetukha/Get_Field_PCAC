#include <complex>

#include "ED_Par.h"
#include "Constants.h"
#include "element_geom.h"

//==========Падающее поле Hinc(x) = e^i(k, x) * (k x eo) / (omega * mu * mu0)====================
void get_hinc(const double* x, const double* e0, const ED_Par& ed_param, std::complex<double>* h_inc)
{
    std::complex<double>  deg, deg1;
    double vvv[3];
    deg = std::complex<double>(1., 0.) / ed_param.omega0 / ed_param.mu_d[0] / Constants::m0;
    vec_prod(ed_param.k_vec, e0, vvv);
    deg1 = deg * exp(std::complex<double>(0., 1.) * scal_prod(ed_param.k_vec, x));
    h_inc[0] = vvv[0] * deg1;
    h_inc[1] = vvv[1] * deg1;
    h_inc[2] = vvv[2] * deg1;
}