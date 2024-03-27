#ifndef _ELEMENT_GEOM_H_
#define _ELEMENT_GEOM_H_

#include <complex>

//====================================Absolute value=============================================

double abs_tmp(const double a);

double abs_tmp(const std::complex<double> a);

//======================================My comparate=============================================

int cmp(const std::complex<double> a, const std::complex<double> b);

int cmp(const std::complex<double> a, const double b);

int cmp(const double a, const std::complex<double> b);

int cmp(const double a, const double b);


//======================================Scalar product===========================================

double scal_prod(const double* vec_1, const double* vec_2);

std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2);

std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2);

std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2);

//=====================================Vector product============================================
void vec_prod(const double* vec_1, const double* vec_2, double* res);

void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res);

void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);

void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);


//====================================Distance btw 2 points======================================
template <typename T>
T dist(const T* vec_1, const T* vec_2)
{
    T res = sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
    return res;
}

//====================================Square of triagnle=========================================
double tr_square(const double* pnt_1, const double* pnt_2, const double* pnt_3);

//======================================Vector length============================================
template <typename T>
double vec_length(const T* vec_1)
{
    return sqrt(abs_tmp(vec_1[0]) * abs_tmp(vec_1[0]) + abs_tmp(vec_1[1]) *
                    abs_tmp(vec_1[1]) + abs_tmp(vec_1[2]) * abs_tmp(vec_1[2]));
}

//=======================================Solid angle=============================================
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x);

//=======================================Cell normal=============================================
void norm_func(const double (&rut0)[4][3], double* norm_res);

//====================Perpendicular to the segment(проекция точки на отрезок)====================
void perp(const double* a, const double* b, const double* x, double* n);


//==================================Nearest point to the edge====================================
//Для статической ячейки
void near_point(const double (&seg)[2][3], const double* x, double* z, double& dist);

void near_point(const double* a, const double* b, const double* x, double* z, double& dist_res);


//=============================Normal to cell edge(нормаль к краю)===============================
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu);

void get_nu(const double* a, const double* b, const double* norm, double* nu);


//=====================================Div Modulo================================================
int div_mod(const int a, const int b);//деление числа a по модулю b(a <= b)

//===================================Cell diameter===============================================
double get_diam(const double (&root_tmp)[4][3]);
double get_diam_triangle(const double (&root_tmp)[3][3]);

//===========================================Center of mass======================================
void get_center_mass(const double* a, const double* b, const double* c, const double* d, double* res);

//===========================================Quadr square========================================
double quadr_square(const double* a, const double* b, const double* c, const double* d);


#endif
