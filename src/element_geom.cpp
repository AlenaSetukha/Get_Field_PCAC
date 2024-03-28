#include <iostream>
#include <cmath>
#include <complex>

#include "constants.h"


//========================================Absolute value=========================================

double abs_tmp(const double a)
{
    return fabs(a);
}

double abs_tmp(const std::complex<double> a)
{
    return sqrt(a.real() * a.real() + a.imag() * a.imag());
}

//========================================My comparate===========================================

int cmp(const std::complex<double> a, const std::complex<double> b)
{
    double a_abs = sqrt(a.real() * a.real() + a.imag() * a.imag());
    double b_abs = sqrt(b.real() * b.real() + b.imag() * b.imag());
    if (a_abs < b_abs)
    {
        return -1;
    } else if (a_abs > b_abs)
    {
        return 1;
    } else return 0;
}

int cmp(const std::complex<double> a, const double b)
{
    double a_abs = sqrt(a.real() * a.real() + a.imag() * a.imag());
    if (a_abs < b)
    {
        return -1;
    } else if (a_abs > b)
    {
        return 1;
    } else return 0;
}

int cmp(const double a, const std::complex<double> b)
{
    double b_abs = sqrt(b.real() * b.real() + b.imag() * b.imag());
    if (a < b_abs)
    {
        return -1;
    } else if (a > b_abs)
    {
        return 1;
    } else return 0;
}

int cmp(const double a, const double b)
{
    if (a < b)
    {
        return -1;
    } else if (a > b)
    {
        return 1;
    } else return 0;
}


//====================================Scalar product=============================================

double scal_prod(const double* vec_1, const double* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const double* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const double* vec_1, const std::complex<double>* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

std::complex<double> scal_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2)
{
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

//====================================Vector product=============================================
void vec_prod(const double* vec_1, const double* vec_2, double* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const std::complex<double>* vec_1, const double* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const double* vec_1, const std::complex<double>* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

void vec_prod(const std::complex<double>* vec_1, const std::complex<double>* vec_2, std::complex<double>* res)
{
    res[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    res[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    res[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
    return;
}

//==================================Distance btw 2 points========================================
template <typename T>
T dist(const T* vec_1, const T* vec_2)
{
    T res = sqrt((vec_1[0] - vec_2[0]) * (vec_1[0] - vec_2[0]) +
        (vec_1[1] - vec_2[1]) * (vec_1[1] - vec_2[1]) +
        (vec_1[2] - vec_2[2]) * (vec_1[2] - vec_2[2]));
    return res;
}
    
//====================================Square of triagnle=========================================
double tr_square(const double* pnt_1, const double* pnt_2, const double* pnt_3)
{
    double ab[3], ac[3];
    for (int i = 0; i < 3; i++)
    {
        ab[i] = pnt_2[i] - pnt_1[i];
        ac[i] = pnt_3[i] - pnt_1[i];
    }
    double n[3];
    vec_prod(ab, ac, n);
    double res = scal_prod(n, n);
    res = sqrt(res) / 2.;
    return res;
}

//======================================Vector length============================================
template <typename T>
double vec_length(const T* vec_1)
{
    return sqrt(abs_tmp(vec_1[0]) * abs_tmp(vec_1[0]) + abs_tmp(vec_1[1]) *
                    abs_tmp(vec_1[1]) + abs_tmp(vec_1[2]) * abs_tmp(vec_1[2]));
}

//=======================================Solid angle=============================================
double solid_angle(const double* x_a, const double* x_b, const double* x_c, const double* x)
{
    double r1[3], r2[3], r3[3], l1, l2, l3, deg, v[3];
    r1[0] = x_a[0] - x[0];
    r1[1] = x_a[1] - x[1];
    r1[2] = x_a[2] - x[2];

    r2[0] = x_b[0] - x[0];
    r2[1] = x_b[1] - x[1];
    r2[2] = x_b[2] - x[2];

    r3[0] = x_c[0] - x[0];
    r3[1] = x_c[1] - x[1];
    r3[2] = x_c[2] - x[2];

    l1 = vec_length(r1);
    l2 = vec_length(r2);
    l3 = vec_length(r3);

    //vec_prod(r2, r3, v);
    vec_prod(r3, r2, v);// = r2 x r3(почему-то перевернуто)

    deg = scal_prod(r1, v) / (l1 * l2 * l3 +
            scal_prod(r1, r2) * l3 + scal_prod(r2, r3) * l1 + scal_prod(r3, r1) * l2);
    return 2 * atan(deg);
}

//=======================================Cell normal=============================================
void norm_func(const double (&rut0)[4][3], double* norm_res)
{
    double ac[3], bd[3];
    for (int i = 0; i < 3; i++)
    {
        ac[i] = rut0[2][i] - rut0[0][i];
        bd[i] = rut0[3][i] - rut0[1][i]; 
    }
    double v[3], l;
    vec_prod(ac, bd, v);
    l = vec_length(v);
    norm_res[0] = v[0] / l;
    norm_res[1] = v[1] / l;
    norm_res[2] = v[2] / l;
}

//======================Perpendicular to the segment(проекция точки на отрезок)==================
void perp(const double* a, const double* b, const double* x, double* n)
{
    double ab[3], ax[3];
    for (int i = 0; i < 3; i++)
    {
        ab[i] = b[i] - a[i];
        ax[i] = x[i] - a[i];
    }
    double sc_pr = scal_prod(ab, ax) / scal_prod(ab, ab);
    for (int i = 0; i < 3; i++)
    {
        n[i] = a[i] + (b[i] - a[i]) * sc_pr;
    }  
    return;
}

//==================================Nearest point to the segment=================================
//статическая функция
void near_point(const double (&seg)[2][3], const double* x, double* z, double& dist_res)
{
    Constants c;
    double a[3], b[3], n[3];
    double dist_ab, dist_an, dist_bn;
    for (int i = 0; i < 3; i++)
    {
        a[i] = seg[0][i];
        b[i] = seg[1][i];
    }
    perp(a, b, x, n);
    dist_ab = dist(a,b);
    dist_an = dist(a,n);
    dist_bn = dist(b,n); 
    if (dist_an + dist_bn - dist_ab < c.machine_zero)
    {
        z[0] = n[0];
        z[1] = n[1];
        z[2] = n[2];
        dist_res = dist(x, n);
    } else {
        if (dist(x, a) < dist(x, b))
        {
            dist_res = dist(x, a);
            z[0] = a[0];
            z[1] = a[1];
            z[2] = a[2];
        } else {
            dist_res = dist(x, b);
            z[0] = b[0];
            z[1] = b[1];
            z[2] = b[2];
        }
    }
    return;
}


void near_point(const double* a, const double* b, const double* x, double* z, double& dist_res)
{
    Constants c;
    double n[3];
    double dist_ab, dist_an, dist_bn;
    perp(a, b, x, n);
    dist_ab = dist(a,b);
    dist_an = dist(a,n);
    dist_bn = dist(b,n); 
    if (dist_an + dist_bn - dist_ab < c.machine_zero)
    {
        z[0] = n[0];
        z[1] = n[1];
        z[2] = n[2];
        dist_res = dist(x, n);
    } else {
        if (dist(x, a) < dist(x, b))
        {
            dist_res = dist(x, a);
            z[0] = a[0];
            z[1] = a[1];
            z[2] = a[2];
        } else {
            dist_res = dist(x, b);
            z[0] = b[0];
            z[1] = b[1];
            z[2] = b[2];
        }
    }
    return;
}


//==============================Normal to cell edge(нормаль к краю)==============================
//rut0 - ячейка с 4 вершинами
//norm - нормаль к ячейке
//i, inext - индексы вершин A,B 
void get_nu(const double (&rut0)[4][3], const double* norm, const int i, const int inext, double* nu)
{
    Constants c;
    double len_nu, diff[3]; 
    for (int j = 0; j < 3; j++)
    {
        diff[j] = rut0[inext][j] - rut0[i][j];
    }
    vec_prod(diff, norm, nu);
    len_nu = vec_length(nu);

    if (len_nu > c.machine_zero)
    {
        for (int j = 0; j < 3; j++)
        {
            nu[j] /= len_nu;
        }
    }
    return;
}

//a, b - концы отрезка, norm - нормаль к ячейке
void get_nu(const double* a, const double* b, const double* norm, double* nu)
{
    Constants c;
    double len_nu, diff[3]; 
    for (int j = 0; j < 3; j++)
    {
        diff[j] = b[j] - a[j];
    }
    vec_prod(diff, norm, nu);
    len_nu = vec_length(nu);

    if (len_nu > c.machine_zero)
    {
        for (int j = 0; j < 3; j++)
        {
            nu[j] /= len_nu;
        }
    }
    return;
}




//==================================Деление по модулю===========================================
int div_mod(const int a, const int b)//деление числа a по модулю b, b <= a
{
    if (a == b)
    {
        return 0;
    } else
    {
        return a;
    }
}

//========================================Cell diameter==========================================
double get_diam(const double (&root_tmp)[4][3])
{
    double diag1[3], diag2[3], res = 0.;
    for (int g = 0; g < 3; g++)
    {
        diag1[g] = root_tmp[2][g] - root_tmp[0][g];
        diag2[g] = root_tmp[3][g] - root_tmp[1][g];
    }
    if (vec_length(diag1) > vec_length(diag2))
    {
        res = vec_length(diag1);
    } else {
        res = vec_length(diag2);
    }


    if (dist(root_tmp[1], root_tmp[0]) > res)
    {
        res = dist(root_tmp[1], root_tmp[0]);
    }
    if (dist(root_tmp[2], root_tmp[1])  > res)
    {
        res = dist(root_tmp[2], root_tmp[1]);
    }
    if (dist(root_tmp[3], root_tmp[2])  > res)
    {
        res = dist(root_tmp[3], root_tmp[2]);
    }
    if (dist(root_tmp[0], root_tmp[3])  > res)
    {
        res = dist(root_tmp[0], root_tmp[3]);
    }
    return res;
}


double get_diam_triangle(const double (&root_tmp)[3][3])
{
    double len1, len2, len3;
    len1 = dist(root_tmp[0], root_tmp[1]);
    len2 = dist(root_tmp[1], root_tmp[2]);
    len3 = dist(root_tmp[2], root_tmp[0]);

    double res = 0.;
    if (len1 > len2)
    {
        res = len1;
    } else {
        res = len2;
    }

    if (len3 > res)
    {
        res = len3;
    }
    return res;
}

//========================================Center of mass=========================================
void get_center_mass(const double* a, const double* b, const double* c, const double* d, double* res)
{
    for (int i = 0; i < 3; i++)
    {
        res[i] = (a[i] + b[i] + c[i] + d[i]) * 0.25;
    }
    return;
}


//========================================Quadr square===========================================
double quadr_square(const double* a, const double* b, const double* c, const double* d)
{
    double ac[3], bd[3], vec[3];
    for (int i = 0; i < 3; i++)
    {
        ac[i] = c[i] - a[i];
        bd[i] = d[i] - b[i];
    }
    vec_prod(ac, bd, vec);
    return vec_length(vec) / 2.;
}

