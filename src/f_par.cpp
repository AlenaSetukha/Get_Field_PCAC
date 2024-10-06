#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "f_par.h"

f_par::f_par()
{
    rs = 0., calc_dist = 0.;
    k = 0.;
    for (int i = 0; i < 3; i++)
    {
        n[i] = 0., vec_cmplx[i] = 0., ort[i] = 0.;
        vec_dbl[i] = 0., e0[i] = 0;
    }
    eps_edge = 0.;
}

f_par::f_par(const double rs_in)
{
    rs = rs_in;

    calc_dist = 0.;
    k = 0.;
    for (int i = 0; i < 3; i++)
    {
        n[i] = 0., vec_cmplx[i] = 0., ort[i] = 0.;
        vec_dbl[i] = 0., e0[i] = 0;
    }
    eps_edge = 0.;
}

f_par::f_par(const double rs_in, const std::complex<double> k_in)
{
    rs = rs_in;
    k = k_in;

    calc_dist = 0.;
    for (int i = 0; i < 3; i++)
    {
        n[i] = 0., vec_cmplx[i] = 0., ort[i] = 0.;
        vec_dbl[i] = 0., e0[i] = 0;
    }
    eps_edge = 0.;
}

f_par::f_par(const double rs_in, const std::complex<double> k_in, const double* n_in)
{
    rs = rs_in;
    k = k_in;
    n[0] = n_in[0];
    n[1] = n_in[1];
    n[2] = n_in[2];


    calc_dist = 0.;
    for (int i = 0; i < 3; i++)
    {
        vec_cmplx[i] = 0., ort[i] = 0.;
        vec_dbl[i] = 0., e0[i] = 0;
    }
    eps_edge = 0.;
}

f_par::f_par(const f_par& obj)
{
    rs = obj.rs;
    calc_dist = obj.calc_dist;
    k = obj.k;
    eps_edge = obj.eps_edge;
    for (int i = 0; i < 3; i++)
    {
        n[i] = obj.n[i], vec_cmplx[i] = obj.vec_cmplx[i];
        vec_dbl[i] = obj.vec_dbl[i], ort[i] = obj.ort[i], e0[i] = obj.e0[i];
    }
}
