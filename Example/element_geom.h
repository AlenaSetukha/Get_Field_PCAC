#ifndef _ELEMENT_GEOM_H_
#define _ELEMENT_GEOM_H_

#include <complex>
//======================================My comparate=============================================

int cmp(const std::complex<double> a, const std::complex<double> b);

int cmp(const std::complex<double> a, const double b);

int cmp(const double a, const std::complex<double> b);

int cmp(const double a, const double b);


//======================================Scalar product===========================================
double scal_prod(const double (&vec_1)[3], const double (&vec_2)[3]);

std::complex<double> scal_prod(const std::complex<double> (&vec_1)[3], const double (&vec_2)[3]);

std::complex<double> scal_prod(const double (&vec_1)[3], const std::complex<double> (&vec_2)[3]);

std::complex<double> scal_prod(const std::complex<double> (&vec_1)[3], const std::complex<double> (&vec_2)[3]);

//=====================================Vector product============================================
void vec_prod(const double (&vec_1)[3], const double (&vec_2)[3], double (&res)[3]);

void vec_prod(const std::complex<double> (&vec_1)[3], const double (&vec_2)[3], std::complex<double> (&res)[3]);

void vec_prod(const double (&vec_1)[3], const std::complex<double> (&vec_2)[3], std::complex<double> (&res)[3]);

void vec_prod(const std::complex<double> (&vec_1)[3], const std::complex<double> (&vec_2)[3], std::complex<double> (&res)[3]);


//====================================Distance btw 2 points======================================
template <typename T>
T dist(const T (&vec_1)[3], const T (&vec_2)[3])
{
    T res = sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
    return res;
}

//====================================Square of triagnle=========================================
double tr_square(const double (&pnt_1)[3], const double (&pnt_2)[3], const double (&pnt_3)[3]);

//======================================Vector length============================================
template <typename T>
double vec_length(const T (&vec_1)[3]);

//=======================================Solid angle=============================================
double solid_angle(const double (&x_a)[3], const double (&x_b)[3], const double (&x_c)[3], const double (&x)[3]);

//=======================================Cell normal=============================================
void norm_func(const double (&rut0)[4][3], double (&norm_res)[3]);

//====================Perpendicular to the segment(проекция точки на отрезок)====================
void perp(const double (&a)[3], const double (&b)[3], const double (&x)[3], double (&n)[3]);


//==================================Nearest point to the edge====================================
//Для статической ячейки
void near_point(const double (&seg)[2][3], const double (&x)[3], double (&z)[3], double& dist_res);

void near_point(const double (&a)[3], const double (&b)[3], const double (&x)[3], double (&z)[3], double& dist_res);


//=============================Normal to cell edge(нормаль к краю)===============================
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu);

void get_nu(const double (&a)[3], const double (&b)[3], const double (&norm)[3], double (&nu)[3]);


//===================================Cell diameter===============================================
double get_diam(const double (&root_tmp)[4][3]);
double get_diam_triangle(const double (&root_tmp)[3][3]);

//===========================================Center of mass======================================
void get_center_mass(const double (&a)[3], const double (&b)[3], const double (&c)[3], const double (&d)[3], double (&res)[3]);

//===========================================Quadr square========================================
double quadr_square(const double (&a)[3], const double (&b)[3], const double (&c)[3], const double (&d)[3]);


#endif
