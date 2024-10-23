#ifndef _ELEMENT_GEOM_H_
#define _ELEMENT_GEOM_H_

#include <complex>
//=====================================================================================
//----------------Библиотека основных геометрических операций--------------------------
//=====================================================================================
/**
 * Описание функций/процедур:
 *      scal_prod - скалярное произведение векторов vec_1, vec_2. Перегружено
 *      vec_prod - векторное произведение векторов vec_1, vec_2. Перегружено
 *      tr_square - площадь треугольника
 *      quadr_square - площадь четурехугольника
 *      solid_angle - вычисление телесного угла(в радианах)
 *      norm_func - вектор нормали к ячейке. Перегружено
 *      perp - проекция точки на отрезок(вектор перпендикуляра)
 *      near_point - точка(основание) перпендикуляра на отрезке.
 *      get_nu - нормаль к краю
 *      get_diam - размер ячейки(максимальная сторона/диагональ). Перегружено
 *      get_center_mass - центр масс ячейки. Перегружено
 *      check_points_match - проверка на совпадение двух точек
 * 
 *      dist - евклидово расстояние между двумя точками/векторами
 *      vec_length - длина вектора
 */



//==Scalar product==
double scal_prod(const double* vec_1, const double* vec_2);
std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2);
std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2);
std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2);

//==Vector product==
void vec_prod(const double* vec_1, const double* vec_2, double* res);
void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res);
void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);
void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res);


//==Square of triagnle==
double tr_square(const double* pnt_1, const double* pnt_2, const double* pnt_3);
//==Quadr square==
double quadr_square(const double* a, const double* b, const double* c, const double* d);


//==Solid angle==
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x);

//==Cell normal==
void norm_func(const double (&rut0)[4][3], double* norm_res);
void norm_func(const double (&rut0)[3][3], double* norm_res);

//==Perpendicular to the segment(проекция точки на отрезок)==
void perp(const double* a, const double* b, const double* x, double* n);


//==Nearest point to the edge==
void near_point(const double (&seg)[2][3], const double* x, double* z, double& dist_res);
void near_point(const double* a, const double* b, const double* x, double* z, double& dist_res);


//==Normal to cell edge(нормаль к краю)==
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu);
void get_nu(const double* a, const double* b, const double* norm, double* nu);


//==Cell diameter==
double get_diam(const double (&root_tmp)[4][3]);
double get_diam(const double (&root_tmp)[3][3]);

//==Center of mass==
void get_center_mass(const double* a, const double* b, const double* c, const double* d, double* res);
void get_center_mass(const double* a, const double* b, const double* c, double* res);
void get_center_mass(const double (&root_tmp)[4][3], double* res);
void get_center_mass(const double (&root_tmp)[3][3], double* res);


//==Check two points are different==
int check_points_match(const double* a, const double* b);




//==Distance btw 2 points(vectors)==
template <typename T>
T dist(const T* vec_1, const T* vec_2)
{
    T res = sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
    return res;
}


//==Vector length==
template <typename T>
double vec_length(const T* vec_1)
{
    return sqrt(std::abs(vec_1[0]) * std::abs(vec_1[0]) +
                std::abs(vec_1[1]) * std::abs(vec_1[1]) +
                std::abs(vec_1[2]) * std::abs(vec_1[2]));
}
#endif