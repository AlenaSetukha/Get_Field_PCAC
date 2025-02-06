#ifndef _GET_HINC_H_
#define _GET_HINC_H_

#include <complex>
#include "ED_Par.h"

//==========Падающее поле Hinc(x) = e^i(k, x) * (k x eo) / (omega * mu * mu0)====================
void get_hinc(const double* x, const double* e0, const ED_Par& ed_param, std::complex<double>* h_inc);

#endif // _GET_HINC_H_