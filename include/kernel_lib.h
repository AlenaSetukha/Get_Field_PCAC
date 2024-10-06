#ifndef _KERNEL_LIB_H_
#define _KERNEL_LIB_H_

#include <complex>
#include "f_par.h"
#include "constants.h"

using namespace Constants;

//=====================================================================================
//------------------Библиотека основных подынтегральных ядер---------------------------
//=====================================================================================
/**
 * Описание функций:
 *      f_simple_pot_L - потенциал простого слоя для ур-я Лапласа
 *      f_grad_simple_pot_L - градиент потенциала простого слоя уравнения Лапласа
 *      f_double_pot_L - потенциал двойного слоя для ур-я Лапласа
 *      f_vector_pot_L - векторный потенциал для уравнения Лапласа
 * 
 *      f_simple_pot_G - потенциал простого слоя для ур-я Гельмгольца
 *      f_double_pot_G - потенциал двойного слоя для ур-я Гельмгольца
 *      f_grad_simple_pot_G - градиент потенциала простого слоя уравнения Гельмгольца
 *      f_grad_simple_pot_G_eps - градиент потенциала простого слоя уравнения Гельмгольца eps
 * 
 * 
 *      f_KFar - ядро для K=rotrot  в дальней зоне
 *      func3 - спецфункция (x - y) / |x - y|^2
 * 
 *      f_aver_simple -ядро осреднения для поверхностной дивергенции(далеко от края)
 *      f_aver_near_edge - ядро осреднения для поверхностной дивергенции(близко к краю)
 *      psi_aver_func - функция осреднения(базовая)
 * 
 *      f_simple_pot_G1 - потенциал простого слоя для ур-я Гельмгольца - 1 
 *      get_einc_ker - падающее поле как ядро Ker = Einc(x) * ort(стороны ячейки) 
 */


void f_simple_pot_L(const double* x, const double* y, const f_par& param, double* ff);

void f_grad_simple_pot_L(const double* x, const double* y, const f_par& param, double* ff);

void f_double_pot_L(const double* x, const double* y, const f_par& param,  double* ff);

template <typename P>
void f_vector_pot_L(const double* x, const double* y, const f_par& param, P* ff)
{
    //F = (1 / 4pi) * a * (y - x) / |x - y|^3
    double r, t, f;
    r = sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) +
            (x[2] - y[2]) * (x[2] - y[2]));

    if (r < Constants::machine_zero)
    {
        for (int i = 0; i < 3; i++)
        {
            ff[i] = static_cast<P>(0);
        }
    } else {
        f = Constants::pi_reverse / r / r / r; //вне эпсилон круга, не сглаживаем функцию, K(x, y)

        double diff[3];
        diff[0] = (y[0] - x[0]) * f;
        diff[1] = (y[1] - x[1]) * f;
        diff[2] = (y[2] - x[2]) * f;

        vec_prod(diff, param.vec_dbl, ff);

        if (r >= Constants::machine_zero && r < param.rs)
        {
            t = r / param.rs;
            for (int i = 0; i < 3; i++)
            {
                ff[i] = ff[i] * (3 * t * t - 2 * t * t * t);//Keps
            }
        }
    }
}



void f_simple_pot_G(const double* x, const double* y, const f_par& param, std::complex<double>* ff);

void f_double_pot_G(const double* x, const double* y, const f_par& param,  std::complex<double>* ff);

void f_grad_simple_pot_G(const double* x, const double* y, const f_par& param, std::complex<double>* ff);

void f_grad_simple_pot_G_eps(const double* x, const double* y, const f_par& param, std::complex<double>* ff);




void f_KFar(const double* x, const double* y, const f_par& param, std::complex<double>* ff);


//--(x - y) / |x - y|^2--
void func3(const double* x, const double* y, const f_par& param, double* ff);



void f_aver_simple(const double* x, const double* y, const f_par& param, double* res);
void f_aver_near_edge(const double* x, const double* y, const f_par& param, double* res);
double psi_aver_func(const double* x, const double* y, const double eps);


void f_simple_pot_G1(const double* x, const double* y, const f_par& param, std::complex<double>* ff);

void get_einc_ker(const double* x, const f_par& param, std::complex<double>* e_inc);

#endif

