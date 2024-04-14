#ifndef _GET_EINC_H_
#define _GET_EINC_H_

#include <complex>

//=======================Падающее поле Einc(x) = e0 * e^i(k, x)==================================
void get_einc(const double (&x)[3], const double (&e0)[3], const double (&k_vec)[3], std::complex<double> (&e_inc)[3]);

#endif