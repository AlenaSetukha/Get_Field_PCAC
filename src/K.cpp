#include <iostream>
#include <complex>

#include "integral_par.h"
#include "f_par.h"
#include "integral_universal_pnt.h"
#include "k0.h"
#include "kernel_lib.h"
#include "integrals_analytic.h"
#include "Num_Par.h"




//-----------------Оператор  K = rotrot в ближней зоне с выделением особенности------------------

void K_rot_rot_Near(const std::complex<double>* j, const double* x,
        const double (&rut0)[4][3], const double* norm,
        const std::complex<double> k, const Num_Par& num_param,
        std::complex<double>* res)
{
    // -curv
    f_par param_seg(num_param.rs_seg * get_diam(rut0) / num_param.n_start_seg, k);
    integral_par  int_parGradF(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);
    k0(j, x, norm, rut0, f_grad_simple_pot_G, param_seg, int_parGradF, res);


    // k^2 с выделением особенности
    f_par param(num_param.rs * get_diam(rut0) / num_param.n_start, k);
    integral_par  int_parF(1, num_param.n_start, num_param.p_max, num_param.eps);
    param.calc_dist = 3. * get_diam(rut0); // когда считать аналитически
    std::complex<double> cur_res;
    integral_simple_pot_G(x, rut0, param, int_parF, cur_res);
    for (int i = 0; i < 3; i++)
    {
        res[i] += cur_res * k * k * j[i];
    }
}
