#ifndef _K0_H_
#define _K0_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "f_par.h"
#include "integral_par.h"
#include "constants.h"
#include "element_geom.h"
#include "integral_universal_seg_pnt.h"

using namespace Constants;


//===============================================================================================
// K0 - сумма 4-х интегралов по сторонам ячейки rut0, специальная функция.
// Параметры:
//      e - вектор, для которого считаем(j)
//      x - точка коллокации
//      norm - нормаль к ячейке
//      rut0 - ячейка(4 вершины)
//      f_0 - функция ядра
//      param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
//      int_param - параметры вычисления интеграла
//      res - результат вычисления суммы четырех инетгралов
//===============================================================================================

template <typename P>
void k0(const P (&e)[3], const double (&x)[3],
        const double (&norm)[3], const double (&rut0)[4][3],
        void (*f_0)(const double(&)[3], const double(&)[3], const f_par&, P*),
        const f_par& param, const integral_par& int_par,
        P* res)
{
    P* ks = new P[int_par.idim];
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < int_par.idim; g++)
    {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    for (int s = 0; s < 4; s++)
    {
        for (int g = 0; g < int_par.idim; g++)
        {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }

        if (s != 3)
        {

            integral_universal_seg_pnt(rut0[s], rut0[s + 1], x, f_0, param, int_par, ks);

            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[s + 1][g] - rut0[s][g];//b-a
            }
        } else {
            integral_universal_seg_pnt(rut0[3], rut0[0], x, f_0, param, int_par, ks);
            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[0][g] - rut0[3][g];//b-a
            }
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);


        if (len_nu > Constants::machine_zero)
        {
            for (int g = 0; g < int_par.idim; g++)
            {
                nu[g] /= len_nu;
            }
            for (int g = 0; g < int_par.idim; g++)
            {
                ks[g] *= -scal_prod(e, nu);//ks[3] ab
                res[g] += ks[g];
            }
        }
    }

    delete[] ks;
}






//===============================================================================================
// K0 - сумма 3-х интегралов по сторонам ячейки rut0, специальная функция.
// Параметры:
//      e - вектор, для которого считаем(j)
//      x - точка коллокации
//      norm - нормаль к ячейке
//      rut0 - ячейка(3 вершины)
//      f_0 - функция ядра
//      param - параметры функции(k, радиус сглаживания, машинный ноль, доп параметры )
//      int_param - параметры вычисления интеграла
//      res - результат вычисления суммы трех инетгралов
//===============================================================================================

template <typename P>
void k0(const P (&e)[3], const double (&x)[3],
        const double (&norm)[3], const double (&rut0)[3][3],
        //void (*f_0)(const double*, const double*, const f_par&, P*),
        void (*f_0)(const double(&)[3], const double(&)[3], const f_par&, P*),
        const f_par& param, const integral_par& int_par,
        P* res)
{
    P* ks = new P[int_par.idim];
    double len_nu;
    double nu[3], diff[3];

    for (int g = 0; g < int_par.idim; g++)
    {
        res[g] = static_cast<P>(0);
        ks[g] = static_cast<P>(0);
    }


    for (int s = 0; s < 3; s++)
    {
        for (int g = 0; g < int_par.idim; g++)
        {
            ks[g] = static_cast<P>(0);
            diff[g] = 0.;
        }

        if (s != 2)
        {
            integral_universal_seg_pnt(rut0[s], rut0[s + 1], x, f_0, param, int_par, ks);

            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[s + 1][g] - rut0[s][g];//b-a
            }
        } else {
            integral_universal_seg_pnt(rut0[3], rut0[0], x, f_0, param, int_par, ks);
            for (int g = 0; g < int_par.idim; g++)
            {
                diff[g] = rut0[0][g] - rut0[3][g];//b-a
            }
        }

        vec_prod(diff, norm, nu);
        len_nu = vec_length(nu);


        if (len_nu > Constants::machine_zero)
        {
            for (int g = 0; g < int_par.idim; g++)
            {
                nu[g] /= len_nu;
            }
            for (int g = 0; g < int_par.idim; g++)
            {
                ks[g] *= -scal_prod(e, nu);//ks[3] ab
                res[g] += ks[g];
            }
        }
    }

    delete[] ks;
}
#endif

