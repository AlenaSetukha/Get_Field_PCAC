#ifndef _K_H_
#define _K_H_

#include <iostream>
#include <complex>

#include "integral_par.h"
#include "f_par.h"
#include "integral_universal_pnt.h"
#include "k0.h"
#include "kernel_lib.h"
#include "integrals_analytic.h"
#include "Num_Par.h"


//===============================================================================================
//------------------------Оператор  K = rotrot в дальней зоне------------------------------------
//===============================================================================================
/**
 * K[sigma[j], j[j]] = surf_int(Kfar). Полностью численно, со сглаживанием.
 * Kfar = (e^(ikr) / 4pi) * (j_vec * f1 + (x - y)(x - y, j_vec) * f2 / r^2),
 *      f1 = -1 / r^3 + (ik) / r^2 + k^2 / r,
 *      f2 = 3 / r^3 - (3ik) / r^2 - k^2 / r.
 * Параметры:
 *      j[3] - ток на ячейке
 *      x[3] - точка расчета
 *      rut0[4/3][3] - ячейка(четырехугольная или треугольная)
 *      k - волновое число подынтегральной функции
 *      num_param - численные параметры(для интегрирования) дальней зоны
 *      res[3] - результирующий вектор
 * 
 * num_param.rs - сглаживание относительно diam текущей ячейки(~0.5-1)
 */


template<size_t CellPoints>
void K_rot_rot_Far(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3],
        const std::complex<double> k,
        const Num_Par& num_param,
        std::complex<double>* res)
{
    // Инициализация параметров
    f_par param(num_param.rs * get_diam(rut0), k);
    integral_par int_parF(3, num_param.n_start, num_param.p_max, num_param.eps);
    param.vec_cmplx[0] = j[0], param.vec_cmplx[1] = j[1], param.vec_cmplx[2] = j[2];

    integral_universal_pnt(x, rut0, f_KFar, param, int_parF, res);
}








//===============================================================================================
//-------------------Оператора K = rotrot в ближней зоне без выделения особенности---------------
//===============================================================================================
/**
 * K[sigma[j], j[j]] = -curv_int(gradx F) + k^2 * surf_int(F), F = eikr/r 
 * Полностью численно, оба интеграла со сглаживанием.
 * Параметры:
 *      j[3] - ток на ячейке
 *      x[3] - точка расчета
 *      rut0[4/3][3] - ячейка(четырехугольная или треугольная)
 *      norm[3] - нормаль к ячейке
 *      num_param - численные параметры(для интегралов)
 *      res[3] - результирующий вектор
 * 
 * rs/rs_seg - радиус сглаживания пов/крив интегралов относительно
 * (diamj / n_start), (~0.5 - 2)
 */

template<size_t CellPoints>
void K_rot_rot(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3], const double* norm,
        const Num_Par& num_param,
        const std::complex<double> k,
        std::complex<double>* res)
{
    // Инициализация параметров
    f_par param(num_param.rs * get_diam(rut0) / num_param.n_start, k);
    f_par param_seg(num_param.rs_seg * get_diam(rut0) / num_param.n_start_seg, k);
    integral_par int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    integral_par  int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);

    //k^2
    std::complex<double> cur_res[1];
    integral_universal_pnt(x, rut0, f_simple_pot_G, param, int_parF, cur_res);

    std::complex<double> tmp = cur_res[0] * param.k * param.k;
    res[0] = tmp * j[0];
    res[1] = tmp * j[1];
    res[2] = tmp * j[2];

    //grad div
    std::complex<double> cur_res3[3];
    k0(j, x, norm, rut0, f_grad_simple_pot_G, param_seg, int_parGradF, cur_res3);

    res[0] += cur_res3[0];
    res[1] += cur_res3[1];
    res[2] += cur_res3[2];
}    









//===============================================================================================
//-----------------Оператор  K = rotrot в ближней зоне с выделением особенности------------------
//===============================================================================================
/**
 * K[sigma[j], tau[j]] = -curv_int(F_eps) + k^2 * j * surf_int(F), F = eikr/r 
 * Интеграл по краю: численно со сглаживанием. 
 * Интеграл по ячейке: выделение особенности + полуаналитически
 * Параметры:
 *      j[3] - ток на ячейке
 *      x[3] - точка расчета
 *      rut0[4][3] - четырехугольная ячейка
 *      norm[3] - нормаль к ячейке
 *      k - волновое число
 *      num_param - численные параметры для интегрирования
 *      res[3] - результирующий вектор
 * 
 * num_param.rs/rs_seg - в долях от  (diamj / n_start) (~0.5-1)
 */

void K_rot_rot_Near(const std::complex<double>* j, const double* x,
        const double (&rut0)[4][3], const double* norm,
        const std::complex<double> k, const Num_Par& num_param,
        std::complex<double>* res);


#endif