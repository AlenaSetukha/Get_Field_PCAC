#ifndef _R_H_
#define _R_H_

#include <iostream>
#include <complex>

#include "integral_par.h"
#include "f_par.h"
#include "kernel_lib.h"
#include "element_geom.h"
#include "Num_Par.h"
#include "constants.h"



//===============================================================================================
//------------------Вычисление оператора R = rot в дальней зоне----------------------------------
//===============================================================================================
/**
 * R[sigma[j], tau[j]] = -j x surf_int(gradF), F = eikr/r
 * Параметры:
 *      j[3] - магнитный ток на ячейке
 *      x[3] - точка расчета
 *      rut0 - ячейка(четыерхугольная или треугольная)
 *      num_param - параметры подынтегральной функции
 *      grid_step - шаг сетки
 *      res - результирующий вектор
 * 
 * num_param.rs - сглаживание для пов. интеграла
 * относительно diam текущей ячейки(~0.5-1)
*/

template<size_t CellPoints>
void R_rot_Far(const std::complex<double>* j, const double* x,
        const double (&rut0)[CellPoints][3],
        const Num_Par& num_param,
        const std::complex<double> k,
        std::complex<double>* res)
{
    std::complex<double> cur_res3[3]{};
    // Инициализация параметров
    f_par param(num_param.rs * get_diam(rut0), k);
    integral_par int_parGradF(3, num_param.n_start, num_param.p_max, num_param.eps);

    integral_universal_pnt(x, rut0, f_grad_simple_pot_G, param, int_parGradF, cur_res3);
    vec_prod(cur_res3, j, res);
}






//===============================================================================================
//-----------------Вычисление оператора R = rot в точках коллокации------------------------------
//===============================================================================================
/**
 * R[sigma[j], tau[j]] (x_i) = -j x surf_int(gradF), F = eikr/r
 * R = 0, если i = j
 * Параметры:
 *      j[3] - магнитный ток на ячейке
 *      x[3] - точка коллокации
 *      rut0 - ячейка(четыерхугольная или треугольная)
 *      num_param - параметры подынтегральной функции
 *      grid_step - шаг сетки
 *      res - результирующий вектор
 * 
 * num_param.rs - сглаживания нет
 */

template<size_t CellPoints>
void R_rot_Colloc(const std::complex<double>* j,
        const double* x,
        const double (&rut0)[CellPoints][3],
        const Num_Par& num_param,
        const std::complex<double> k,
        std::complex<double>* res)
{
    double y[3];
    get_center_mass(rut0, y);
    if (dist(x, y) < Constants::machine_zero) // i = j
    {
        res[0] = 0., res[1] = 0., res[2] = 0.;
    } else {
        std::complex<double> cur_res3[3]{};
        // Инициализация параметров
        f_par param(Constants::machine_zero, k);
        integral_par int_parGradF(3, num_param.n_start, num_param.p_max, num_param.eps);

        integral_universal_pnt(x, rut0, f_grad_simple_pot_G, param, int_parGradF, cur_res3);
        vec_prod(cur_res3, j, res);
    }
}
#endif