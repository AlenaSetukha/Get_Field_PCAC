#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "constants.h"
#include "element_geom.h"
#include "f_par.h"

static inline double sqr(const double x) { return x * x; }

//=====================================================================================
//-----------------------Потенциал простого слоя для ур-я Лапласа----------------------
//=====================================================================================
void f_simple_pot_L(const double *x, const double *y, const f_par &param,
                    double *ff) {
  // F = 1 / (4pi * |x - y|)
  double r;
  r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));
  if (r < Constants::machine_zero) {
    ff[0] = 0.0;
  } else {
    ff[0] = 1. / (4 * M_PI) / r;
    if (r >= Constants::machine_zero && r < param.rs) {
      double t = r / param.rs;
      ff[0] = ff[0] * (3 * t * t - 2 * t * t * t);
    }
  }
}

//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Лапласа-------------------
//=====================================================================================

void f_grad_simple_pot_L(const double *x, const double *y, const f_par &param,
                         double *ff) {
  // F = (y - x) / (4pi * |x - y|^3)
  double r;
  r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));
  if (r < Constants::machine_zero) {
    for (int i = 0; i < 3; i++) {
      ff[i] = 0.;
    }
  } else {
    for (int i = 0; i < 3; i++) {
      ff[i] = (y[i] - x[i]) * 1. / (4 * M_PI) / r / r / r;
    }
    if (r < param.rs) {
      double t;
      t = r / param.rs;
      for (int i = 0; i < 3; i++) {
        ff[i] *= (3 * t * t - 2 * t * t * t);
      }
    }
  }
}

//=====================================================================================
//---------------------Потенциал двойного слоя для ур-я Лапласа------------------------
//=====================================================================================
void f_double_pot_L(const double *x, const double *y, const f_par &param,
                    double *ff) {
  // F = (1 / 4pi) * (n(y) x (x - y) / |x - y|^3)
  double r;
  r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));

  if (r < Constants::machine_zero) {
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

    ff[0] = 1. / (4 * M_PI) * scal_prod(norm, diff) / r / r / r; // K(x,y)
    if (r < param.rs) {
      double t = r / param.rs;
      ff[0] = ff[0] * (3 * t * t - 2 * t * t * t); // Keps(x,y)
    }
  }
}

//=====================================================================================
//---------------------Потенциал простого слоя для ур-я Гельмгольца--------------------
//=====================================================================================
void f_simple_pot_G(const double *x, const double *y, const f_par &param,
                    std::complex<double> *ff) {
  // F = (1 / 4pi) * (e^ik|x - y| / |x - y|)
  double r = sqrt(sqr(x[0] - y[0]) + sqr(x[1] - y[1]) + sqr(x[2] - y[2]));
  ff[0] = (1. / (4 * M_PI)) * exp(std::complex<double>(0., r) * param.k);
  double t = r / param.rs;
  ff[0] *= (t < 1) * (t * (3 - 2 * t) / param.rs) + (t >= 1) / r;
}

//=====================================================================================
//-----------------------Потенциал двойного слоя для ур-я Гельмгольца------------------
//=====================================================================================
void f_double_pot_G(const double *x, const double *y, const f_par &param,
                    std::complex<double> *ff) {
  // F = (1 / 4pi) * (n(y) x (x - y) / |x - y|^3) * (e^ik|x - y| - ike^ik|x -
  // y|)
  double r =
      sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));

  if (r < Constants::machine_zero) {
    ff[0] = std::complex<double>(0., 0.);
  } else {
    double diff[3];
    diff[0] = x[0] - y[0];
    diff[1] = x[1] - y[1];
    diff[2] = x[2] - y[2];

    double norm[3];
    norm[0] = param.n[0];
    norm[1] = param.n[1];
    norm[2] = param.n[2];

    ff[0] = 1. / (4 * M_PI) * scal_prod(norm, diff) / (r * r * r);

    std::complex<double> ikr = std::complex<double>(0., 1.) * param.k * r;
    std::complex<double> e_ikr = exp(ikr);
    std::complex<double> ike_ikr =
        std::complex<double>(0., 1.) * param.k * exp(ikr);

    ff[0] = ff[0] * (e_ikr - ike_ikr);
    if (r < param.rs) {
      double t = r / param.rs;
      ff[0] = ff[0] * t * t * (3 - 2 * t);
    }
  }
}

//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Гельмгольца---------------
//=====================================================================================
void f_grad_simple_pot_G(const double *x, const double *y, const f_par &param,
                         std::complex<double> *ff) {
  // F = (1 / 4pi) * (ikr - 1) * e^ikr * (x - y) / r^3
  double r;
  std::complex<double> f;

  r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
           (x[2] - y[2]) * (x[2] - y[2]));

  if (r < Constants::machine_zero) {
    for (int i = 0; i < 3; i++) {
      ff[i] = std::complex<double>(0., 0.);
    }
  } else {
    f = 1. / (4 * M_PI) / (r * r * r);
    std::complex<double> ikr = std::complex<double>(0., 1.) * param.k * r;
    f *= exp(ikr) * (ikr - std::complex<double>(1., 0.));
    for (int i = 0; i < 3; i++) {
      ff[i] = (x[i] - y[i]) * f;
    }
    if (r < param.rs) {
      double t = r / param.rs;
      for (int i = 0; i < 3; i++) {
        ff[i] *= (t * t * (3 - 2 * t));
      }
    }
  }
}

//=====================================================================================
//---------------Градиент потенциала простого слоя уравнения Гельмгольца eps-----------
//=====================================================================================
void f_grad_simple_pot_G_eps(const double *x, const double *y,
                             const f_par &param, std::complex<double> *ff) {
  double r, r0, tau[3];

  tau[0] = param.vec_dbl[0];
  tau[1] = param.vec_dbl[1];
  tau[2] = param.vec_dbl[2];
  std::complex<double> f;

  r0 = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));
  r = sqrt(r0 * r0 + param.rs * param.rs);

  if (r0 < Constants::machine_zero || r < Constants::machine_zero) {
    for (int i = 0; i < 3; i++) {
      ff[i] = std::complex<double>(0., 0);
    }
  } else {
    double diff[3];
    diff[0] = (x[0] - y[0]) / r0;
    diff[1] = (x[1] - y[1]) / r0;
    diff[2] = (x[2] - y[2]) / r0;

    std::complex<double> f1 =
        (-1.0 / r / r / r) + (std::complex<double>(0., 1.) * param.k / r / r);
    std::complex<double> f2 =
        (3.0 / r / r / r) -
        (3.0 * std::complex<double>(0., 1.) * param.k / r / r) -
        (param.k * param.k / r);
    f2 *= scal_prod(diff, tau);

    std::complex<double> eikr =
        exp(std::complex<double>(0., 1.) * param.k * r0) * 1. /
        (4 * M_PI); // попробовать e^ikr

    ff[0] = (diff[0] * f2 + tau[0] * f1) * eikr;
    ff[1] = (diff[1] * f2 + tau[1] * f1) * eikr;
    ff[2] = (diff[2] * f2 + tau[2] * f1) * eikr;
  }
}

//=====================================================================================
//---------------------Ядро для оператора K(rotrot) в дальней зоне---------------------
//=====================================================================================
void f_KFar(const double *x, const double *y, const f_par &param,
            std::complex<double> *ff) {
  /**
   * Kfar = (e^(ikr) / 4pi) * (j_vec * f1 + (x - y)(x - y, j_vec) * f2 / r^2)
   * f1 = -1 / r^3 + (ik) / r^2 + k^2 / r
   * f2 = 3 / r^3 - (3ik) / r^2 - k^2 / r
   */

  double r = dist(x, y);

  if (r < Constants::machine_zero) {
    for (int i = 0; i < 3; i++) {
      ff[i] = std::complex<double>(0., 0);
    }
  } else {
    std::complex<double> j[3];
    j[0] = param.vec_cmplx[0], j[1] = param.vec_cmplx[1],
    j[2] = param.vec_cmplx[2];
    double diff[3];
    for (int i = 0; i < 3; i++) {
      diff[i] = (x[i] - y[i]) / r;
    }
    std::complex<double> f1 =
        (-1.0 / r / r / r) + (std::complex<double>(0., 1.) * param.k / r / r +
                              param.k * param.k / r);
    std::complex<double> f2 =
        (3.0 / r / r / r) -
        (3.0 * std::complex<double>(0., 1.) * param.k / r / r) -
        (param.k * param.k / r);
    f2 *= scal_prod(diff, j);
    std::complex<double> eikr =
        exp(std::complex<double>(0., 1.) * param.k * r) * 1. / (4 * M_PI);
    for (int i = 0; i < 3; i++) {
      ff[0] = (diff[i] * f2 + j[i] * f1) * eikr;
    }

    if (r < param.rs) {
      double t = r / param.rs;
      for (int i = 0; i < 3; i++) {
        ff[i] *= (t * t * (3 - 2 * t));
      }
    }
  }
}

//=====================================================================================
//-------------------Ядро для интеграла от (x - y) / |x - y|^2-------------------------
//=====================================================================================
void func3(const double *x, const double *y, const f_par &param, double *ff) {
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