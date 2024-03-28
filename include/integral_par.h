#ifndef _INTEGRAL_PAR_H_
#define _INTEGRAL_PAR_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

class integral_par
{
public:
    int idim;    // размерность подынтегральной функции
    int n_start; // стартовое разбиение ячейки/отрезка
    int p_max;   // предельное число разбиений(итераций)- 2^p_max
    double eps;  // точность вычисления интеграла(имеет смысл при p_max != 1)

    integral_par(const int idim_in, const int n_start_in, const int p_max_in, const double eps_in);

    integral_par(const integral_par& obj);
};
#endif

