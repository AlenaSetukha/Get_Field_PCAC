#ifndef _INTEGRAL_UNIVERSAL_SEG_PNT_H_
#define _INTEGRAL_UNIVERSAL_SEG_PNT_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "Kernel_Par.h"
#include "Integral_Par.h"
#include "element_geom.h"

//=====================================Description================================================
// Расчет интеграла от функции F(xk,y) по отрезку (в точке). Дифференцирование по y
// Аргументы:
//     a - начало отрезка интегрирования
//     b - конец отрезка интегрирования
//     x - точка, в которой считаем инетграл(точка коллокации)
//     f_0 - функция ядра F(xk, y)
//     param - параметры функции ядра(k, радиус сглаживания, машинный ноль, доп параметры )
//     int_param - параметры интегрирования(разбиение, размерность, и тд)
//     res - результат вычисления интеграла
//================================================================================================

template <typename P>
void integral_universal_seg_pnt(const double* a, const double* b, const double* x,
        void (*f_0)(const double*, const double*, const Kernel_Par&, P*),
        const Kernel_Par& param, const Integral_Par& int_param,
        P* res)
{
    int p_n;
    int n = int_param.n_start;
    double d[3], y[3], dl, delta = 0.;

    P* ff = new P[int_param.idim];
    P* res_prev = new P[int_param.idim];

    for (int g = 0; g < int_param.idim; g++) {
        ff[g] = static_cast<P>(0);
        res[g] = static_cast<P>(0);
        res_prev[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.p_max; p_n++) {
        for (int j = 0; j < 3; j++) {
            d[j] = (b[j] - a[j]) / n;
        }
        dl = vec_length(d);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                y[j] = a[j] + (i + 0.5) * d[j];
            }
            f_0(x, y, param, ff);

            for (int g = 0; g < int_param.idim; g++) {
                res[g] += ff[g] * static_cast<P>(dl);
            }
        }

        delta = 0.;
        for (int g = 0; g < int_param.idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.eps * int_param.eps && p_n != 0) {
            break;
        }

        n = n * 2;
        for (int g = 0; g < int_param.idim; g++) {
            res_prev[g] = res[g];
            res[g] = static_cast<P>(0);
        }
    }

    if (p_n == int_param.p_max) {
        for (int g = 0; g < int_param.idim; g++) {
            res[g] = res_prev[g];
        }
    }

    delete[] res_prev;
    delete[] ff;
}
#endif // _INTEGRAL_UNIVERSAL_SEG_PNT_H_
