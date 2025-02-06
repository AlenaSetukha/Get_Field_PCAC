#ifndef _INTEGRALS_ANALYTIC_H_
#define _INTEGRALS_ANALYTIC_H_

#include <array>

//==============================================================================================
//--------------------------------Integral 1/|x-y|----------------------------------------------
//==============================================================================================
//Брать эту функцию, или численную - решать вне.
double integral1Divr(const double (&rut0)[4][3], const double* x);
double integral1Divr(const double (&rut0)[3][3], const double *x);



//==============================================================================================
//--------------------------------Integral d/dn(1/|x-y|)----------------------------------------
//==============================================================================================
//Четырехугольная ячейка, точка не лежит на ячейке
void integral_ddn_1Divr(const double (&rut0)[4][3], const double* x, double& res);



//==============================================================================================
//--------------------------------Integral (x-y)/(|x-y|^2)--------------------------------------
//==============================================================================================
//Ниже предполагается четырехугольная ячейка, точка лежит на ячейке
void integral_xmyDivr2(const double (&rut0)[4][3], const double* x, double* res);



//==============================================================================================
//--------------------------------Integral nu/|x-y|(curvilinear)--------------------------------
//==============================================================================================
//Четырехугольная ячейка
void integralnu1Divr(const double (&rut0)[4][3], const double* x, const double eps, double* res);





//==============================================================================================
//--------------Полуаналитический поверхностный интеграл по ячейке от функции-------------------
//-----------------------F(x-y) = e^ik|x - y| / 4 * pi * |x - y|--------------------------------
//==============================================================================================
#include "Kernel_Par.h"
#include "Integral_Par.h"
#include "integral_universal_pnt.h"
#include "kernel_lib.h"
#include "element_geom.h"


//=====================================================================================
//------------------Функция, которая возвращает полуаналитический----------------------
//-------------------поверхностный интеграл по ячейке от функции-----------------------
//---------------------F(x-y) = e^ik|x - y| / 4 * pi * |x - y|-------------------------
//=====================================================================================
template<size_t CellPoints>
void integral_simple_pot_G(const double *x, const double (&rut0)[CellPoints][3],
                           const Kernel_Par &param, const Integral_Par &int_param,
                           std::complex<double> &ff)
{
    double y[3];
    get_center_mass(rut0, y);
    std::complex<double> tmp[1];
    ff = std::complex<double>(0., 0.);

    if (dist(x, y) < param.calc_dist) {
        // близко, считаем аналитически, разбивая на 2 интеграла
        ff = static_cast<std::complex<double>>(integral1Divr(rut0, x)) * 1. /
             (4 * M_PI);
        if (dist(x, y) < 1e-5) {
            ff += cell_square(rut0) *
                  std::complex<double>(0., 1.) * param.k / (4 * M_PI);
        } else {
            integral_universal_pnt(x, rut0, f_simple_pot_G1, param, int_param, tmp);
            ff += tmp[0];
        }
    } else {
        // далеко, считаем интеграл численно
        integral_universal_pnt(x, rut0, f_simple_pot_G, param, int_param, tmp);
        ff = tmp[0];
    }
}

#endif // _INTEGRALS_ANALYTIC_H_