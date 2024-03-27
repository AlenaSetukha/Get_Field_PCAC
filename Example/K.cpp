#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>


#include "integral_par.h"
#include "f_par.h"
#include "integral_universal_pnt.h"
#include "k0.h"
#include "kernel_lib.h"


void K_rot_rot(const std::complex<double> (&j)[3], const double* x,
        const double (&rut0)[4][3], const double* norm,
        const integral_par& integral_par_f_simple, const integral_par& integral_par_f_grad_simple,
        const f_par& param, const f_par& param_seg,
        std::complex<double>* res)
{
    //k^2
    std::complex<double> cur_res[1];
    integral_universal_pnt(x, rut0, f_simple_pot_G, param, integral_par_f_simple, cur_res);

    std::complex<double> tmp = cur_res[0] * param.k * param.k;
    res[0] = tmp * j[0];
    res[1] = tmp * j[1];
    res[2] = tmp * j[2];

    //grad div
    std::complex<double> cur_res3[3];
    k0(j, x, norm, rut0, f_grad_simple_pot_G, param_seg, integral_par_f_grad_simple, cur_res3);

    res[0] += cur_res3[0];
    res[1] += cur_res3[1];
    res[2] += cur_res3[2];

    return;
}




void K_rot_rot(const std::complex<double> (&j)[3], const double* x,
        const double (&rut0)[3][3], const double* norm,
        const integral_par& integral_par_f_simple, const integral_par& integral_par_f_grad_simple,
        const f_par& param, const f_par& param_seg,
        std::complex<double>* res)
{
    //k^2
    std::complex<double> cur_res[1];
    integral_universal_pnt(x, rut0, f_simple_pot_G, param, integral_par_f_simple, cur_res);

    std::complex<double> tmp = cur_res[0] * param.k * param.k;
    res[0] = tmp * j[0];
    res[1] = tmp * j[1];
    res[2] = tmp * j[2];

    //grad div
    std::complex<double> cur_res3[3];
    k0(j, x, norm, rut0, f_grad_simple_pot_G, param_seg, integral_par_f_grad_simple, cur_res3);

    res[0] += cur_res3[0];
    res[1] += cur_res3[1];
    res[2] += cur_res3[2];

    return;
}
