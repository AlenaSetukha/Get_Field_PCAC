#ifndef _GET_EINC_H_
#define _GET_EINC_H_

#include <complex>

//=======================Падающее поле Einc(x) = e0 * e^i(k, x)==================================
void get_einc(const double* x, const double* e0, const double* k_vec, std::complex<double>* e_inc);

#endif // _GET_EINC_H_