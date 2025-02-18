#ifndef _K0_H_
#define _K0_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "Kernel_Par.h"
#include "Integral_Par.h"
#include "element_geom.h"
#include "integral_universal_seg_pnt.h"
#include "Constants.h"



//===============================================================================================
//-----------------Оператор  K0 - криволинейный интеграл-----------------------------------------
//===============================================================================================
/**
 * K0 - криволинейный интеграл по контуру ячейки от функции вида
 *              K0 = int_curv { (j, nu) * F(x-y) dy}
 *                  где F(x - y) - некоторое ядро
 * Параметры:
 *      j - вектор, для которого считаем интеграл
 *      x - точка коллокации
 *      rut0 - ячейка(4 или 3 вершины)
 *      f_0 - функция ядра
 *      param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
 *      int_param - параметры вычисления интеграла
 *      res - результат вычисления суммы четырех инетгралов
 */
//===============================================================================================
template <typename P, size_t CellPoints>
void K0(const P* j, const double* x,
    const double (&rut0)[CellPoints][3],
    void (*f_0)(const double*, const double*, const Kernel_Par&, P*),
    const Kernel_Par& param, const Integral_Par& int_par,
    P* res)
{
    P* ks = new P[int_par.idim];
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < int_par.idim; g++) {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    double norm[3];
    norm_func(rut0, norm);

    for (int s = 0; s < (int)CellPoints; s++) {
        for (int g = 0; g < int_par.idim; g++) {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }
        if (s != (int)CellPoints - 1) {
            integral_universal_seg_pnt(rut0[s], rut0[s + 1], x, f_0, param, int_par, ks);

            for (int g = 0; g < int_par.idim; g++) {
                diff[g] = rut0[s + 1][g] - rut0[s][g];//b-a
            }
        } else {
            integral_universal_seg_pnt(rut0[(int)CellPoints - 1], rut0[0], x, f_0, param, int_par, ks);
            for (int g = 0; g < int_par.idim; g++) {
                diff[g] = rut0[0][g] - rut0[(int)CellPoints - 1][g];//b-a
            }
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);


        if (len_nu > Constants::machine_zero) {
            for (int g = 0; g < int_par.idim; g++) {
                nu[g] /= len_nu;
            }
            for (int g = 0; g < int_par.idim; g++) {
                ks[g] *= scal_prod(j, nu);//ks[3] ab
                res[g] += ks[g];
            }
        }
    }
    delete[] ks;
}
#endif // _K0_H_

