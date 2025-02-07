#ifndef _K_H_
#define _K_H_

#include <iostream>
#include <complex>

#include "Integral_Par.h"
#include "Kernel_Par.h"
#include "integral_universal_pnt.h"
#include "K0.h"
#include "kernel_lib.h"
#include "integrals_analytic.h"
#include "Num_Par.h"


#include "element_geom.h"

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
 * num_param.rs - радиус сглаживания интеграла
 * относительно размера ячеки второго уровня (~0.5 - 2) 
 */

template<size_t CellPoints>
void K_rot_rot_Far(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3],
        const Num_Par& num_param, const std::complex<double> k, 
        std::complex<double>* res)
{
    // Инициализация параметров
    Kernel_Par param_KFar(num_param.rs * get_diam(rut0) / num_param.n_start, k);
    Integral_Par int_parF(3, num_param.n_start, num_param.p_max, num_param.eps);

    param_KFar.vec_cmplx[0] = j[0], param_KFar.vec_cmplx[1] = j[1], param_KFar.vec_cmplx[2] = j[2];

    integral_universal_pnt(x, rut0, f_KFar, param_KFar, int_parF, res);
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
 *      num_param - численные параметры(для интегралов)
 *      res[3] - результирующий вектор
 * 
 * rs/rs_seg - радиус сглаживания пов/крив интегралов
 * относительно размера ячеки второго уровня (~0.5 - 2) 
 */

template<size_t CellPoints>
void K_rot_rot(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3],
        const Num_Par& num_param, const std::complex<double> k,
        std::complex<double>* res)
{
    // Инициализация параметров
    Kernel_Par param_F(num_param.rs * get_diam(rut0) / num_param.n_start, k);
    Kernel_Par param_GradF(num_param.rs_seg * get_diam(rut0) / num_param.n_start_seg, k);
    Integral_Par int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    Integral_Par  int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);

    //k^2
    std::complex<double> cur_res[1];
    integral_universal_pnt(x, rut0, f_simple_pot_G, param_F, int_parF, cur_res);


    std::complex<double> tmp = cur_res[0] * param_F.k * param_F.k;
    res[0] = tmp * j[0];
    res[1] = tmp * j[1];
    res[2] = tmp * j[2];

    //grad div
    std::complex<double> cur_res3[3];
    K0(j, x, rut0, f_grad_simple_pot_G, param_GradF, int_parGradF, cur_res3);
    res[0] -= cur_res3[0];
    res[1] -= cur_res3[1];
    res[2] -= cur_res3[2];
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
 *      k - волновое число
 *      num_param - численные параметры для интегрирования
 *      res[3] - результирующий вектор
 * 
 * num_param.rs/rs_seg - радиус сглаживания пов/крив интегралов
 * относительно размера ячеки второго уровня (~0.5 - 2)
 */
template<size_t CellPoints>
void K_rot_rot_Near(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3],
        const Num_Par& num_param, const std::complex<double> k,
        std::complex<double>* res)
{
    // k^2 с выделением особенности
    Kernel_Par param(num_param.rs * get_diam(rut0) / num_param.n_start, k);
    Integral_Par  int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    param.calc_dist = 3. * get_diam(rut0); // когда считать аналитически
    
    std::complex<double> cur_res;
    integral_simple_pot_G(x, rut0, param, int_parF, cur_res);
    for (int i = 0; i < 3; i++) {
        res[i] = cur_res * k * k * j[i];
    }

    // -curv
    std::complex<double> res3[3];
    Kernel_Par param_seg(num_param.rs_seg * get_diam(rut0) / num_param.n_start_seg, k);
    Integral_Par  int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    K0(j, x, rut0, f_grad_simple_pot_G, param_seg, int_parGradF, res3);
    for (int i = 0; i < 3; i++) {
        res[i] -= res3[i];
    }
}

#endif // _K_H_
