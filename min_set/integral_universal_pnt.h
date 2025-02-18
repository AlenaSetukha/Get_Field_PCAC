#ifndef _INTEGRAL_UNIVERSAL_PNT_H_
#define _INTEGRAL_UNIVERSAL_PNT_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "Kernel_Par.h"
#include "Integral_Par.h"
#include "element_geom.h"


//=====================================Description================================================
// Расчет поверхностного интеграла по четырехугольной ячейке от функции
// F(xk, y) в точке xk. Дифференцирование по y.
// Аргументы:
//     x - точка коллокации(точка, в которой считается инетграл)
//     rut0 - ячейка(4 вершины по 3 координаты)
//     f_0 - функция ядра F(xk, y)
//     param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
//     int_param - параметры интегрирования(разбиение, размерность, и тд)
//     res - результат вычисления интеграла
//================================================================================================
template<typename P>
void integral_universal_pnt(const double* x, const double (&rut0)[4][3],
        void (*f_0)(const double*, const double*, const Kernel_Par&, P*),
        const Kernel_Par& param, const Integral_Par& int_param, P* res)
{
    double p, q, p1, q1, s, a[3], b[3], a1[3], a2[3], a3[3], a4[3], m1[3], m2[3], rc[3], rn[3];
    double delta = 0.;
    int n = int_param.n_start;
    int p_n;

    P* ff = new P[int_param.idim]();
    P* res_prev = new P[int_param.idim]();
    for (int g = 0; g < int_param.idim; g++) {
        res[g] = static_cast<P>(0);
    }


    for (p_n = 0; p_n < int_param.p_max; p_n++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p = (double)i / (double)n;
                q = (double)j / (double)n;
                p1 = ((double)i + 1.) / (double)n;
                q1 = ((double)j + 1.) / (double)n;

                for (int k = 0; k < 3; k++) {
                    a[k] = q * rut0[1][k] + (1 - q) * rut0[0][k];
                    b[k] = q * rut0[2][k] + (1 - q) * rut0[3][k];
                    a1[k] = p * b[k] + (1 - p) * a[k];
                    a4[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    a[k] = q1 * rut0[1][k] + (1 - q1) * rut0[0][k];
                    b[k] = q1 * rut0[2][k] + (1 - q1) * rut0[3][k];
                    a2[k] = p * b[k] + (1 - p) * a[k];
                    a3[k] = p1 * b[k] + (1 - p1) * a[k];
                }

                for (int k = 0; k < 3; k++) {
                    rc[k] = (a1[k] + a2[k] + a3[k] + a4[k]) / 4.0;
                    m1[k] = ((a2[k] + a3[k]) - (a1[k] + a4[k])) / 2.0;
                    m2[k] = ((a3[k] + a4[k]) - (a1[k] + a2[k])) / 2.0;
                }
                vec_prod(m1, m2, rn);
                s = vec_length(rn);

                f_0(x, rc, param, ff);

                for (int g = 0; g < int_param.idim; g++) {
                    res[g] += ff[g] * static_cast<P>(s);
                }
            }
        }

        delta = 0.;
        for (int g = 0; g < int_param.idim; g++) {
            delta += std::abs(res[g] - res_prev[g]) * std::abs(res[g] - res_prev[g]);
        }

        if (delta <  int_param.eps * int_param.eps && p_n != 0) {
            break;
        }

        n = n * 2; //шаг уменьшается в 2 раза
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










//=====================================Description================================================
// Расчет поверхностного интеграла по треугольной ячейке от функции
// F(xk, y) в точке xk. Дифференцирование по y.
// Аргументы:
//     x - точка коллокации(точка, в которой считается инетграл)
//     rut0 - ячейка(3 вершины по 3 координаты)
//     f_0 - функция ядра F(xk, y)
//     param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
//     int_param - параметры интегрирования(разбиение, размерность, и тд)
//     res - результат вычисления интеграла
//================================================================================================
template<typename P>
void integral_universal_pnt(const double* x, const double (&rut0)[3][3],
        void (*f_0)(const double*, const double*, const Kernel_Par&, P*),
        const Kernel_Par& param, const Integral_Par& int_param, 
        P* res)
{
    double s, rc[3], A_jP[3], B_j[3], C_jP[3], A_jM[3], C_jM[3];
    double p_vec[3], q_vec[3], r_vec[3];
    double delta = 0.;
    int n =  int_param.n_start, p_n;


    P* ff = new P[int_param.idim];
    P* res_prev = new P[int_param.idim];
    for (int g = 0; g < int_param.idim; g++) {
        ff[g] = static_cast<P>(0);
        res[g] = static_cast<P>(0);
        res_prev[g] = static_cast<P>(0);
    }

    s = tr_square(rut0[0], rut0[1], rut0[2]); // площадь всего треугольника


    for (p_n = 0; p_n < int_param.p_max; p_n++) {
        for (int k = 0; k < 3; k++) {
            p_vec[k] = (rut0[1][k] - rut0[0][k]) / n; // p
            q_vec[k] = (rut0[2][k] - rut0[0][k]) / n; // q
            r_vec[k] = (rut0[2][k] - rut0[1][k]) / n; // r
        }
        for (int i = 0; i < n; i++) { // m
            for (int j = 0; j <= n - 1 - i; j++) { // k
                for (int k = 0; k < 3; k++) {
                    B_j[k] = rut0[0][k] + (i + 1) * p_vec[k] + j * q_vec[k];
                    A_jM[k] = B_j[k] - p_vec[k];
                    C_jM[k] = A_jM[k] + q_vec[k];
                    rc[k] = (B_j[k] + A_jM[k] + C_jM[k]) / 3.0;
                }

                f_0(x, rc, param, ff);
                for (int g = 0; g < int_param.idim; g++) {
                    res[g] += ff[g];
                }

                if (j < n - 1 - i) {
                    for (int k = 0; k < 3; k++) {
                        A_jP[k] = B_j[k] + r_vec[k];
                        C_jP[k] = B_j[k] + q_vec[k];
                        rc[k] = (B_j[k] + A_jP[k] + C_jP[k]) / 3.0;
                    }

                    f_0(x, rc, param, ff);
                    for (int g = 0; g < int_param.idim; g++) {
                        res[g] += ff[g];
                    }
                }
            }
        }

	    delta = 0.;
        for (int g = 0; g < int_param.idim; g++) {
            res[g] = res[g] * s / (double)n / (double)n;
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

    delete[] ff;
    delete[] res_prev;
}
#endif // _INTEGRAL_UNIVERSAL_PNT_H_