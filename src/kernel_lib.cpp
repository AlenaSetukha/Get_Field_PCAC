#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>

#include "constants.h"
#include "element_geom.h"
#include "f_par.h"


//=====================================================================================
//-----------------------Потенциал простого слоя для ур-я Лапласа----------------------
//=====================================================================================
void f_simple_pot_L(const double* x, const double* y, const f_par& param, double* ff)
{
//F = 1 / (4pi * |x - y|)
    Constants c;
    double r, t;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));
    if (r < c.machine_zero)
    {
        ff[0] = 0.0;
    } else {
        ff[0] = c.pi_reverse / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)
        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);//Keps(x, y)
        }
    }
    return;
}


//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Лапласа-------------------
//=====================================================================================

void f_grad_simple_pot_L(const double* x, const double* y, const f_par& param, double* ff)
{
//F = (y - x) / (4pi * |x - y|^3)
    Constants c;
    double r, t;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));
    if (r < c.machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = 0.;
        }
    } else {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = (y[i] - x[i]) * c.pi_reverse / r / r / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)
        }
        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            for (int i = 0; i < 3; i++)
            {
                ff[i] *= (3 * t * t - 2 * t * t * t);//Keps(x, y)
            }
        }
    }
    return;
}


//=====================================================================================
//---------------------Потенциал двойного слоя для ур-я Лапласа------------------------
//=====================================================================================
void f_double_pot_L(const double* x, const double* y, const f_par& param,  double* ff)
{
//F = (1 / 4pi) * (n(y) x (x - y) / |x - y|^3)
    Constants c;
    double r, t;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));

    if (r < c.machine_zero)
    {
        ff[0] = 0.;
    } else {
        double diff[3];
        diff[0] = x[0] - y[0];
        diff[1] = x[1] - y[1];
        diff[2] = x[2] - y[2];

        double norm[3];
        norm[0] = param.n[0];
        norm[1] = param.n[1];
        norm[2] = param.n[2];


        ff[0] = c.pi_reverse * scal_prod(norm, diff) / r / r / r;//K(x,y)
        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);//Keps(x,y)
        }
    }
    return;
}



//=====================================================================================
//---------------------Векторный потенциал для уравнения Лапласа-----------------------
//=====================================================================================
template <typename P>
void f_vector_pot_L(const double* x, const double* y, const f_par& param, P* ff)
{
//F = (1 / 4pi) * a * (y - x) / |x - y|^3

    Constants c;
    double r, t, f;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));

    if (r < c.machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = static_cast<P>(0);
        }
    } else {
        f = c.pi_reverse / r / r / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)

        double diff[3];
        diff[0] = (y[0] - x[0]) * f;
        diff[1] = (y[1] - x[1]) * f;
        diff[2] = (y[2] - x[2]) * f;

        vec_prod(diff, param.a, ff);

        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            for (int i = 0; i < 3; i++)
            {
                ff[i] = ff[i] * (3 * t * t - 2 * t * t * t);//Keps
            }
        }
    }
    return;
}














//=====================================================================================
//---------------------Потенциал простого слоя для ур-я Гельмгольца--------------------
//=====================================================================================
void f_simple_pot_G(const double* x, const double* y, const f_par& param, std::complex<double>* ff)
{
//F = (1 / 4pi) * (e^ik|x - y| / |x - y|)

    Constants c;
    double r;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));


    if (r < c.machine_zero)
    {
        ff[0] = std::complex<double>(0., 0.);
    } else {
        std::complex<double> ikr = c.i_complex * param.k * r;//ikr
        ff[0] = c.pi_reverse * exp(ikr) / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)
        if (r >= c.machine_zero && r < param.rs)
        {
            double t = r / param.rs;
            ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);//Keps(x, y)
        }
    }

    return;
}


//=====================================================================================
//-----------------------Потенциал двойного слоя для ур-я Гельмгольца------------------
//=====================================================================================
void f_double_pot_G(const double* x, const double* y, const f_par& param,  std::complex<double>* ff)
{
//F = (1 / 4pi) * (n(y) x (x - y) / |x - y|^3) * (e^ik|x - y| - ike^ik|x - y|)
    Constants c;
    double r, t;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));

    if (r < c.machine_zero)
    {
        ff[0] = std::complex<double>(0., 0.);
    } else {
        double diff[3];
        diff[0] = x[0] - y[0];
        diff[1] = x[1] - y[1];
        diff[2] = x[2] - y[2];

        double norm[3];
        norm[0] = param.n[0];
        norm[2] = param.n[1];
        norm[3] = param.n[2];

        ff[0] = c.pi_reverse * scal_prod(norm, diff) / r / r / r;

        std::complex<double> ikr = c.i_complex * param.k * r;//ikr
        std::complex<double> e_ikr = exp(ikr);
        std::complex<double> ike_ikr = c.i_complex * param.k * exp(ikr);

        ff[0] = ff[0] * (e_ikr - ike_ikr);//K(x, y)
        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);//Keps(x,y)
        }
    }
    return;
}


//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Гельмгольца---------------
//=====================================================================================
void f_grad_simple_pot_G(const double* x, const double* y, const f_par&param, std::complex<double>* ff)
{
//F = (1 / 4pi) * (ikr - 1) * e^ikr * (x - y) / r^3

    Constants c;
    double r, t;
    std::complex<double> f;

    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));

    if (r < c.machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = std::complex<double>(0., 0.);
        }
    } else {
        f = c.pi_reverse / r / r / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)
        std::complex<double> ikr = c.i_complex * param.k * r;//ikr
        f *= exp(ikr) * (ikr - std::complex<double>(1., 0.));
        for (int i = 0; i < 3; i++)
        {
            ff[i] = (x[i] - y[i]) * f;
        }
        if (r >= c.machine_zero && r < param.rs)
        {
            t = r / param.rs;
            for (int i = 0; i < 3; i++)
            {
                ff[i] = ff[i] * (3 * t * t - 2 * t * t * t);//Keps
            }
        }
    }
    return;
}


//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Гельмгольца eps-----------
//=====================================================================================
void f_grad_simple_pot_G_eps(const double* x, const double* y, const f_par& param, std::complex<double>* ff)
{
    Constants c;
    double r, r0, tau[3];

    tau[0] = param.a[0];
    tau[1] = param.a[1];
    tau[2] = param.a[2];
    std::complex<double> f;

    r0 = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));
    r = sqrt(r0 * r0 + param.rs * param.rs);

    if (r0 < c.machine_zero || r < c.machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = std::complex<double>(0., 0);
        }
    } else {
        double diff[3];
        diff[0] = (x[0] - y[0]) / r0;
        diff[1] = (x[1] - y[1]) / r0;
        diff[2] = (x[2] - y[2]) / r0;

        std::complex<double> f1 = (-1.0 / r / r / r) + (c.i_complex * param.k / r / r);
        std::complex<double> f2 = (3.0 / r / r / r) - (3.0 * c.i_complex * param.k / r / r) - (param.k * param.k / r);
        f2 *= scal_prod(diff, tau);

        std::complex<double> eikr = exp(c.i_complex * param.k * r0) * c.pi_reverse;//попробовать e^ikr

        ff[0] = (diff[0] * f2 + tau[0] * f1) * eikr;
        ff[1] = (diff[1] * f2 + tau[1] * f1) * eikr;
        ff[2] = (diff[2] * f2 + tau[2] * f1) * eikr;
    }
    return;
}







//=====================================================================================
//-------------------Ядро для интеграла от (x - y) / |x - y|^2-------------------------
//=====================================================================================
void func3(const double* x, const double* y, const f_par& param, double* ff)
{
    Constants c;
    double diff[3];
    for (int i = 0; i < 3; i++)
    {
        diff[i] = x[i] - y[i];
    }
    if (vec_length(diff) < c.machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = 0.;
        }
    } else {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = diff[i] / vec_length(diff) / vec_length(diff);
        }
    }
    return; 
}


//=====================================================================================
//---------------Ядро осреднения для поверхностной дивергенции(далеко от края)---------
//=====================================================================================
void f_aver_simple(const double* x, const double* y, const f_par& param, double* res)
{
    double eps_edge = param.eps_edge;
    Constants c;
    double reps = dist(x, y) / eps_edge;
    double deg = (-6.0 + 2. * reps * reps) * exp(-reps * reps) / c.pi / eps_edge / eps_edge / eps_edge / eps_edge; 
    for (int i = 0; i < 3; i++)
    {
        res[i] = (x[i] - y[i]) * deg; 
    }
    return;
}

//=====================================================================================
//---------------Ядро осреднения для поверхностной дивергенции(близко к краю)----------
//=====================================================================================
void f_aver_near_edge(const double* x, const double* y, const f_par& param, double* res)
{
    double eps_edge = param.eps_edge;
    Constants c;
    double reps = dist(x, y) / eps_edge;//эпс из формулы, не сглаживание
    double deg =  (-20.0 + 8. * reps * reps) * exp(-reps * reps) / c.pi / eps_edge / eps_edge / eps_edge / eps_edge; 
    for (int i = 0; i < 3; i++)
    {
        res[i] = (x[i] - y[i]) * deg; 
    }
    return;
}



//=====================================================================================
//---------------------------Функция осреднения(базовая)-------------------------------
//=====================================================================================
//psi(r = |x - y|) = e^(r^2 / eps^2) / (pi * eps^2)
double psi_aver_func(const double* x, const double* y, const double eps)
{
    Constants c;
    double res, deg;

    if (dist(x, y) > eps)
    {
        return 0.;
    }
    deg = dist(x, y) * dist(x, y) / eps / eps;
std::cout << "deg: " << deg << std::endl;
    res = exp(deg) / c.pi / eps / eps;
std::cout << "res: " << res << std::endl;
    return res;
}








//=====================================================================================
//-----------------Потенциал простого слоя для ур-я Гельмгольца - 1--------------------
//=====================================================================================
void f_simple_pot_G1(const double* x, const double* y, const f_par& param, std::complex<double>* ff)
{
//F = (1 / 4pi) * ((e^ik|x - y| - 1) / |x - y|)

    Constants c;
    double r;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));


    if (r < c.machine_zero)
    {
        ff[0] = std::complex<double>(0., 0.);
    } else {
        std::complex<double> ikr = c.i_complex * param.k * r;//ikr
        ff[0] = c.pi_reverse * (exp(ikr) - 1.) / r;//вне эпсилон круга, не сглаживаем функцию, K(x, y)
        if (r >= c.machine_zero && r < param.rs)
        {
            double t = r / param.rs;
            ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);//Keps(x, y)
        }
    }

    return;
}







//===============================================================================================
//-------------------Функция расчета ядра Ker = Einc(x) * ort(стороны ячейки)--------------------
//===============================================================================================
void get_einc_ker(const double* x, const f_par& param, std::complex<double>* e_inc)
{
    Constants c;
    std::complex<double> deg = pow(c.e, c.i_complex * scal_prod(param.a, x));         // a - k_vec, почти всегда double
    e_inc[0] = scal_prod(param.e0, param.ort) * deg;                                  // e0 - double       
}