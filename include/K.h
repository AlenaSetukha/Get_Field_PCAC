#ifndef _K_H_
#define _K_H_

#include <complex>
#include "integral_par.h"
#include "f_par.h"


//===============================================================================================
//-------------------------Вычисление оператора K = rotrot---------------------------------------
//===============================================================================================
// K[sigma[j], tau[j]] = -curv_int() + k^2 * surf_int()m, F = eikr/r
// Параметры:
//      j - подаваемый ток
//      x - точка расчета
//      rut0 - ячейка(4 или 3)
//      norm - нормаль к ячейке
//      integral_par_f_simple/integral_par_f_grad_simple - параметры ингетрирования
//      param/param_seg - параметры подынтегральных функций
//      res - результирующий вектор


void K_rot_rot(const std::complex<double> (&j)[3], const double (&x)[3],
        const double (&rut0)[4][3], const double (&norm)[3],
        const integral_par& integral_par_f_simple, const integral_par& integral_par_f_grad_simple,
        const f_par& param, const f_par& param_seg,
        std::complex<double>* res);

void K_rot_rot(const std::complex<double> (&j)[3], const double (&x)[3],
        const double (&rut0)[3][3], const double (&norm)[3],
        const integral_par& integral_par_f_simple, const integral_par& integral_par_f_grad_simple,
        const f_par& param, const f_par& param_seg,
        std::complex<double>* res);
        
#endif
