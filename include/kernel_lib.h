#ifndef _KERNEL_LIB_H_
#define _KERNEL_LIB_H_

#include <complex>
#include "f_par.h"

//-----------------------Потенциал простого слоя для ур-я Лапласа----------------------
void f_simple_pot_L(const double (&x)[3], const double (&y)[3], const f_par& param, double* ff);

//---------------Градиент потенциала простого слоя уравнения Лапласа-------------------
void f_grad_simple_pot_L(const double (&x)[3], const double (&y)[3], const f_par& param, double* ff);

//---------------------Потенциал двойного слоя для ур-я Лапласа------------------------
void f_double_pot_L(const double (&x)[3], const double (&y)[3], const f_par& param,  double* ff);

//---------------------Векторный потенциал для уравнения Лапласа-----------------------
template <typename P>
void f_vector_pot_L(const double (&x)[3], const double (&y)[3], const f_par& param, P* ff);




//---------------------Потенциал простого слоя для ур-я Гельмгольца--------------------
void f_simple_pot_G(const double (&x)[3], const double (&y)[3], const f_par& param, std::complex<double>* ff);

//-----------------------Потенциал двойного слоя для ур-я Гельмгольца------------------
void f_double_pot_G(const double (&x)[3], const double (&y)[3], const f_par& param,  std::complex<double>* ff);

//---------------Градиент потенциала простого слоя уравнения Гельмгольца---------------
void f_grad_simple_pot_G(const double (&x)[3], const double (&y)[3], const f_par& param, std::complex<double>* ff);

//---------------Градиент потенциала простого слоя уравнения Гельмгольца eps-----------
void f_grad_simple_pot_G_eps(const double (&x)[3], const double (&y)[3], const f_par& param, std::complex<double>* ff);




//-------------------------------(x - y) / |x - y|^2-----------------------------------
void func3(const double (&x)[3], const double (&y)[3], const f_par& param, double* ff);



//---------------Ядро осреднения для поверхностной дивергенции(далеко от края)---------
void f_aver_simple(const double (&x)[3], const double (&y)[3], const f_par& param, double* res);

//---------------Ядро осреднения для поверхностной дивергенции(близко к краю)----------
void f_aver_near_edge(const double (&x)[3], const double (&y)[3], const f_par& param, double* res);



//---------------------------Функция осреднения(базовая)-------------------------------
double psi_aver_func(const double (&x)[3], const double (&y)[3], const double eps);



//-----------------Потенциал простого слоя для ур-я Гельмгольца - 1--------------------
void f_simple_pot_G1(const double (&x)[3], const double (&y)[3], const f_par& param, std::complex<double>* ff);



//------------Падающее поле как ядро Ker = Einc(x) * ort(стороны ячейки)---------------
void get_einc_ker(const double (&x)[3], const f_par& param, std::complex<double>* e_inc);

#endif

