#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "constants.h"
#include "element_geom.h"
#include "f_par.h"

static inline double sqr(const double x) { return x * x; }

/**
 * Функция сглаживания: teta(r / rs) = 3t^2 - 2t^3
 */ 
//=====================================================================================
//-----------------------Потенциал простого слоя для ур-я Лапласа----------------------
//=====================================================================================
void f_simple_pot_L(const double *x, const double *y,
        const f_par &param, double *ff)
{ 
    // F = 1 / (4pi * |x - y|)
    double r = dist(x, y);
    ff[0] = 1. / (4. * M_PI);
    double t = r / param.rs;
    ff[0] *= (t < 1) * (t * (3 - 2 * t) / param.rs) + (t >= 1) / r;
}



//=====================================================================================
//------------------Градиент потенциала простого слоя ур-я Лапласа---------------------
//=====================================================================================
void f_grad_simple_pot_L(const double *x, const double *y,
        const f_par &param, double *ff)
{
    // F = (y - x) / (4pi * |x - y|^3)
    double r = dist(x, y);
    double t = r / param.rs;
    for (int i = 0; i < 3; i++)
    {
        ff[i] = (y[i] - x[i]) / (4. * M_PI);
        ff[i] *= (t < 1) * (3 - 2 * t) / (r * sqr(param.rs)) + (t >= 1) / (sqr(r) * r);
    }
}



//=====================================================================================
//---------------------Потенциал двойного слоя для ур-я Лапласа------------------------
//=====================================================================================
void f_double_pot_L(const double *x, const double *y, const f_par &param, double *ff)
{
    // F = (1 / 4pi) * (n(y) * (x - y) / |x - y|^3)
    double r = dist(x, y);
    double diff[3], norm[3];
    diff[0] = x[0] - y[0], diff[1] = x[1] - y[1], diff[2] = x[2] - y[2];
    norm[0] = param.n[0], norm[1] = param.n[1], norm[2] = param.n[2];
    
    ff[0] = scal_prod(norm, diff) / (4. * M_PI);
    double t = r / param.rs;
    ff[0] *= (t < 1) * (3 - 2 * t) / (r * sqr(param.rs)) + (t >= 1) / (sqr(r) * r);
}










//=====================================================================================
//---------------------Потенциал простого слоя для ур-я Гельмгольца--------------------
//=====================================================================================
void f_simple_pot_G(const double *x, const double *y, const f_par &param,
        std::complex<double> *ff)
{
    // F = (1 / 4pi) * (e^ik|x - y| / |x - y|)
    double r = dist(x, y);
    ff[0] = exp(std::complex<double>(0., r) * param.k) / (4. * M_PI);
    double t = r / param.rs;
    ff[0] *= (t < 1) * (t * (3 - 2 * t) / param.rs) + (t >= 1) / r;
}



//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Гельмгольца---------------
//=====================================================================================
void f_grad_simple_pot_G(const double *x, const double *y, const f_par &param,
        std::complex<double> *ff)
{
    // F = (1 / 4pi) * (ikr - 1) * e^ikr * (x - y) / r^3
    double r = dist(x, y);
    double diff[3];
    diff[0] = x[0] - y[0], diff[1] = x[1] - y[1], diff[2] = x[2] - y[2];

    double t = r / param.rs;
    for (int i = 0; i < 3; i++)
    {
        ff[i] = diff[i] * exp(std::complex<double>(0., r) * param.k) / (4. * M_PI);
        ff[i] *= (std::complex<double>(0., r) * param.k - 1.);
        ff[i] *= (t < 1) * (3 - 2 * t) / (r * sqr(param.rs)) + (t >= 1) / (sqr(r) * r);
    }
}



//=====================================================================================
//-----------------------Потенциал двойного слоя для ур-я Гельмгольца------------------
//=====================================================================================
void f_double_pot_G(const double *x, const double *y, const f_par &param,
        std::complex<double> *ff)
{
    // F = (1 / 4pi) * (n(y) * (y - x) / |x - y|^3) * (e^ikr(ikr - 1))
    double r = dist(x, y);
    double diff[3], norm[3];
    diff[0] = y[0] - x[0], diff[1] = y[1] - x[1], diff[2] = y[2] - x[2];
    norm[0] = param.n[0], norm[1] = param.n[1], norm[2] = param.n[2];

    ff[0] = scal_prod(norm, diff) * exp(std::complex<double>(0., r) * param.k) / (4. * M_PI);
    ff[0] *= (std::complex<double>(0., r) * param.k - 1.);
    double t = r / param.rs;
    ff[0] *= (t < 1) * (3 - 2 * t) / (r * sqr(param.rs)) + (t >= 1) / (sqr(r) * r);
}











//=====================================================================================
//---------------------Ядро для оператора K(rotrot) в дальней зоне---------------------
//=====================================================================================
void f_KFar(const double *x, const double *y, const f_par &param,
      std::complex<double> *ff)
{
/**
* Kfar = (e^(ikr) / 4pi) * (j_vec * f1 + (x - y)(x - y, j_vec) * f2 / r^2)
* f1 = -1 / r^3 + (ik) / r^2 + k^2 / r
* f2 = 3 / r^3 - (3ik) / r^2 - k^2 / r
*/

    double r = dist(x, y);
    if (r < Constants::machine_zero) {
        for (int i = 0; i < 3; i++) {
          ff[i] = std::complex<double>(0., 0.);
        }
    } else {
        std::complex<double> j[3];
        j[0] = param.vec_cmplx[0], j[1] = param.vec_cmplx[1], j[2] = param.vec_cmplx[2];
        double diff[3];
        for (int i = 0; i < 3; i++) {
          diff[i] = (x[i] - y[i]) / r;
        }
        std::complex<double> f1 = (-1.0 / (sqr(r) * r)) + (std::complex<double>(0., 1.) *
                param.k / sqr(r)) + param.k * param.k / r;
        std::complex<double> f2 = (3.0 / (sqr(r) * r)) - (3.0 * std::complex<double>(0., 1.) *
                param.k / sqr(r)) - (param.k * param.k / r);
        f2 *= scal_prod(diff, j);

        double t = r / param.rs;
        for (int i = 0; i < 3; i++) {
            ff[i] = (diff[i] * f2 + j[i] * f1) * exp(std::complex<double>(0., r) * param.k) / (4. * M_PI);
            ff[i] *= (t < 1) * sqr(t) * (3 - 2 * t) + (t >= 1);
        }
    }
}

//=====================================================================================
//-------------------Ядро для интеграла от (x - y) / |x - y|^2-------------------------
//=====================================================================================
void func3(const double *x, const double *y, const f_par &param, double *ff)
{
    double diff[3];
    for (int i = 0; i < 3; i++) {
        diff[i] = x[i] - y[i];
    }
    if (vec_length(diff) < Constants::machine_zero) {
        for (int i = 0; i < 3; i++) {
            ff[i] = 0.;
        }
    } else {
        for (int i = 0; i < 3; i++) {
            ff[i] = diff[i] / vec_length(diff) / vec_length(diff);
        }
    }
}

//=====================================================================================
//---------------Ядро осреднения для поверхностной дивергенции(далеко от края)---------
//=====================================================================================
void f_aver_simple(const double *x, const double *y, const f_par &param,
                   double *res) {
  double eps_edge = param.eps_edge;
  double reps = dist(x, y) / eps_edge;
  double deg = (-6.0 + 2. * reps * reps) * exp(-reps * reps) / M_PI / eps_edge /
               eps_edge / eps_edge / eps_edge;
  for (int i = 0; i < 3; i++) {
    res[i] = (x[i] - y[i]) * deg;
  }
}

//=====================================================================================
//---------------Ядро осреднения для поверхностной дивергенции(близко к краю)----------
//=====================================================================================
void f_aver_near_edge(const double *x, const double *y, const f_par &param,
                      double *res) {
  double eps_edge = param.eps_edge;
  double reps = dist(x, y) / eps_edge; // эпс из формулы, не сглаживание
  double deg = (-20.0 + 8. * reps * reps) * exp(-reps * reps) / M_PI /
               eps_edge / eps_edge / eps_edge / eps_edge;
  for (int i = 0; i < 3; i++) {
    res[i] = (x[i] - y[i]) * deg;
  }
}

//=====================================================================================
//---------------------------Функция осреднения(базовая)-------------------------------
//=====================================================================================
// psi(r = |x - y|) = e^(r^2 / eps^2) / (pi * eps^2)
double psi_aver_func(const double *x, const double *y, const double eps) {
  double res, deg;

  if (dist(x, y) > eps) {
    return 0.;
  }
  deg = dist(x, y) * dist(x, y) / eps / eps;
  return exp(deg) / M_PI / eps / eps;
}

//=====================================================================================
//-----------------Потенциал простого слоя для ур-я Гельмгольца - 1--------------------
//=====================================================================================
void f_simple_pot_G1(const double *x, const double *y, const f_par &param,
                     std::complex<double> *ff) {
  // F = (1 / 4pi) * ((e^ik|x - y| - 1) / |x - y|)

  double r;
  r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));

  if (r < Constants::machine_zero) {
    ff[0] = std::complex<double>(0., 0.);
  } else {
    std::complex<double> ikr = std::complex<double>(0., 1.) * param.k * r;
    ff[0] = 1. / (4 * M_PI) * (exp(ikr) - 1.) / r;
    if (r >= Constants::machine_zero && r < param.rs) {
      double t = r / param.rs;
      ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);
    }
  }
}

//===============================================================================================
//-------------------Функция расчета ядра Ker = Einc(x) * ort(стороны ячейки)--------------------
//===============================================================================================
void get_einc_ker(const double *x, const f_par &param,
                  std::complex<double> *e_inc) {
  std::complex<double> deg =
      std::exp(std::complex<double>(0., 1.) *
               scal_prod(param.vec_dbl, x)); // a - k_vec, почти всегда double
  e_inc[0] = scal_prod(param.e0, param.ort) * deg; // e0 - double
}