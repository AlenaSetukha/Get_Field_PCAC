#ifndef _F_PAR_H_
#define _F_PAR_H_

#include <complex>

class f_par
{
public:
    double rs;              //радиус сглаживания функции ядра
    std::complex<double> k; // волновое число, если есть
    double n[3];            //вектор внешней нормали(если есть/нужен)
    double a[3];            // вектор(a, tau,...)
    double e0[3];
    double ort[3];

    double eps_edge;
    double calc_dist;       //критерий аналитического вычисления, = 3h

    f_par(const double rs_in);

    f_par(const double rs_in, const std::complex<double> k_in);

    f_par(const double rs_in, const std::complex<double> k_in, const double (&n_in)[3]);
};
#endif
