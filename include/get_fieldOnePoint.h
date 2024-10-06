#ifndef _GET_FIELD_ONE_POINT_H_
#define _GET_FILED_ONE_POINT_H_

#include <complex>
#include <vector>

#include "integral_par.h"
#include "f_par.h"
#include "K.h"
#include "element_geom.h"
#include "get_einc.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "constants.h"
#include "R.h"

//===================================================================================
//------------Электрическое поле идеального проводника в дальней зоне----------------
//===================================================================================
/**
 * Расчет полного эл. поля в одной точке с помощью оператора K_far[sigma, j](x):
 *        E_tot(x) = sum_j K_far[sigma_j, j_j](x) + E_inc(x)
 * Параметры:
 *      x - точка расчета
 *      e0[3] - вектор поляризации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      j_vec[num_frm][3] - вектор токов на ячейках в точках коллокаций
 *      cells[num_frm][3/4][3] - список ячеек(четырехугольных/треугольных)
 *      num_par - численные параметры(для интегрирования) дальней зоны
 *      grid_step - шаг сетки
 *      field_E[3] - результат, поле в точке 
 * 
 * num_param.rs - радиус сглаживания (интеграл по ячейке)
 * относительно diam каждой ячейки (0.5 - 1)
 */

template<size_t CellPoints>
void get_fieldOnePoint_idealFar(const double* x, const double* e0,
        const std::complex<double> k0, const double* k_vec,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::vector<double[CellPoints][3]> &cells,
        const Num_Par& num_par,
        std::complex<double>* field_E)
{
    int num_frm = cells.size();
    field_E[0] = 0., field_E[1] = 0., field_E[2] = 0.;
    // Вычисление e_inc, h_inc
    get_einc(x, e0, k_vec, field_E);

    // Вычисление поля 
    std::complex<double> cur_res3[3];
    for (int j = 0; j < num_frm; j++)
    {
        K_rot_rot_Far(j_vec[j], x, cells[j], k0, num_par, cur_res3);
        field_E[0] += cur_res3[0];
        field_E[1] += cur_res3[1];
        field_E[2] += cur_res3[2];
    }
}




//===================================================================================
//---------Электрическое поле идеального проводника в дальней зоне, уточнение--------
//===================================================================================
/**
 * Полное эл. поле в одной точке с путочнением:
 *        E_tot(x) = sum_j K[sigma_j, j_j](x) + E_inc(x)
 * K = K_far или K
 * Параметры:
 *      x - точка расчета
 *      e0[3] - вектор поляризации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      j_vec[num_frm][3] - вектор токов на ячейках в точках коллокаций
 *      cells[num_frm][3/4][3] - список ячеек(четырехугольных/треугольных)
 *      num_parFar - численные параметры(для интегрирования) дальней зоны
 *      num_parNear - параметры ближней зоны без выделения особенности
 *      field_E[3] - результат, поле в точке 
 * 
 * num_parFar.rs - радиус сглаживания (интеграл по ячейке)
 * относительно diam каждой ячейки(0.5 - 1)
 * 
 * num_parNear.rs, rs_seg - радиус сглаживания (интеграл по ячейке/отрезку)
 * относительно (diamj / n_start) (~0.5-2)
 */

template<size_t CellPoints>
void get_fieldOnePoint_idealFarClar(const double* x, const double* e0,
        const std::complex<double> k0, const double* k_vec,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const Num_Par& num_parFar, const Num_Par& num_parNear,
        const double grid_step,
        std::complex<double>* field_E)
{
    int num_frm = cells.size();
    field_E[0] = 0., field_E[1] = 0., field_E[2] = 0.;

    // Вычисление e_inc, h_inc
    get_einc(x, e0, k_vec, field_E);

    // Вычисление поля 
    std::complex<double> cur_res3[3];
    double y[3];

    for (int j = 0; j < num_frm; j++)
    {
        get_center_mass(cells[j], y);
        if (dist(x, y) < 3. * grid_step)
        {
            K_rot_rot(j_vec[j], x, cells[j], norm[j], num_parNear, k0, cur_res3);
        } else {
            K_rot_rot_Far(j_vec[j], x, cells[j], k0, num_parFar, cur_res3);
        }
        field_E[0] += cur_res3[0];
        field_E[1] += cur_res3[1];
        field_E[2] += cur_res3[2];
    }
}





//===================================================================================
//------------Электрическое поле идеального проводника в точках коллокаций-----------
//===================================================================================
/**
 * Полное эл. поле в точке коллокации:
 *        E_tot(x) = sum_j K[sigma_j, j_j](x) + E_inc(x)
 * K = K_far или K_near со сглаживанием
 * Параметры:
 *      x - точка коллокации
 *      e0[3] - вектор поляризации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      j_vec[num_frm][3] - вектор токов на ячейках в точках коллокаций
 *      cells[num_frm][3/4][3] - список ячеек(четырехугольных/треугольных)
 *      num_parFar - численные параметры(для интегрирования) дальней зоны
 *      num_parNear - параметры ближней зоны без выделения особенности
 *      field_E[3] - результат, поле в точке 
 * 
 * num_parFar.rs - радиус сглаживания (интеграл по ячейке)
 * относительно diam каждой ячейки(0.5 - 1)
 * 
 * num_parNear.rs, rs_seg - радиус сглаживания (интеграл по ячейке/отрезку)
 * относительно (diamj / n_start) (~0.5-2)
 */

template<size_t CellPoints>
void get_fieldOnePoint_idealColloc(const double* x, const double* e0,
        const std::complex<double> k0, const double* k_vec,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const Num_Par& num_parFar, const Num_Par& num_parNear,
        const double grid_step,
        std::complex<double>* field_E)
{
    int num_frm = cells.size();
    field_E[0] = 0., field_E[1] = 0., field_E[2] = 0.;

    // Вычисление e_inc, h_inc
    get_einc(x, e0, k_vec, field_E);

    // Вычисление поля 
    std::complex<double> cur_res3[3];
    double y[3];

    for (int j = 0; j < num_frm; j++)
    {
        get_center_mass(cells[j], y);
        if (dist(x, y) < 3. * grid_step)
        {
            K_rot_rot_Near(j_vec[j], x, cells[j], norm[j], k0, num_parNear, cur_res3);
        } else {
            K_rot_rot_Far(j_vec[j], x, cells[j], k0, num_parFar, cur_res3);
        }
        field_E[0] += cur_res3[0];
        field_E[1] += cur_res3[1];
        field_E[2] += cur_res3[2];
    }
}























//===================================================================================
//--------Электрическое поле в одной точке для одиночного диэлектрика----------------
//===================================================================================
/**
 * В зависимости от положения точки: 
 *      get_fieldDielectricA1 - внешняя область диэлектрика,
 *      get_fieldDielectricA2 - внутренняя область диэлектрика.
 * Расчет полного поля в одной точке с помощью операторов K[sigma, j](x), R[sigma, j](x) 
 * для одного облучаемого диэлектрика:
 *      E_tot(x) = (i/(omega eps1)) K1[Sigma, j_E](x) - R1[Sigma, j_M] + E_inc(x) -  внешняя область
 *      E_tot(x) = (i/(omega eps2)) K2[Sigma, j_E](x) - R2[Sigma, j_M] -  внутренняя область
 * Здесь K_rot_rot - ближняя зона, без выделения сообенности,
 *      R_rot - дальняя зона(обычное сглаживание)
 * 
 * 
 * Параметры:
 *      x - точка расчета
 *      e0 - вектор поляризации
 *      num_frm - число ячеек разбиения
 *      cell_list[num_frm] - список четырехугольных ячеек
 *      j_E/j_M [num_frm][3] - вектор эл/магн токов на ячейках в точках коллокаций
 *      ed_param - электродинамические параметры задачи
 *      num_param - численные параметры задачи 
 *      field_E[3] - результат, поле в точке 
 * 
 * 
 * 
 * num_param.rs - сглаживание относительно каждой ячейки(~0.5-1)
 */
template<size_t CellPoints>
void get_fieldDielectricA1(const double* x, const double* e0,
        const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const std::vector<std::complex<double>[3]> &j_E,
        const std::vector<std::complex<double>[3]> &j_M,
        const ED_Par& ed_param, const Num_Par& num_param,
        std::complex<double>* field_E)
{
    field_E[0] = 0., field_E[1] = 0., field_E[2] = 0.;
    
    //=Вычисление поля=
    std::complex<double> sum_K[3]{}, sum_R[3]{}, e_inc[3]{}, cur_res3[3]{};
    std::complex<double> deg = Constants::i_complex /
            (ed_param.omega0 * ed_param.eps_d[0] * Constants::eps0);

    int num_frm = cells.size();
    for (int i = 0; i < num_frm; i++)
    {
        K_rot_rot(j_E[i], x, cells[i], norm[i], num_param, ed_param.k[0], cur_res3);
        sum_K[0] += cur_res3[0];
        sum_K[1] += cur_res3[1];
        sum_K[2] += cur_res3[2];

        R_rot_Far(j_M[i], x, cells[i], num_param, ed_param.k[0], cur_res3);
        sum_R[0] += cur_res3[0];
        sum_R[1] += cur_res3[1];
        sum_R[2] += cur_res3[2];
    }
    get_einc(x, e0, ed_param.k_vec, e_inc);
    field_E[0] = deg * sum_K[0] - sum_R[0] + e_inc[0];
    field_E[1] = deg * sum_K[1] - sum_R[1] + e_inc[1];
    field_E[2] = deg * sum_K[2] - sum_R[2] + e_inc[2];
}





template<size_t CellPoints>
void get_fieldDielectricA2(const double* x,
        const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const std::vector<std::complex<double>[3]> &j_E,
        const std::vector<std::complex<double>[3]> &j_M,
        const ED_Par& ed_param, const Num_Par& num_param,
        std::complex<double>* field_E)
{
    field_E[0] = 0., field_E[1] = 0., field_E[2] = 0.;
    //=Вычисление поля=
    std::complex<double> sum_K[3]{}, sum_R[3]{}, cur_res3[3]{}, j_E_min[3], j_M_min[3];
    std::complex<double> deg = Constants::i_complex /
            (ed_param.omega0 * ed_param.eps_d[1] * Constants::eps0);

    int num_frm = cells.size();
    for (int i = 0; i < num_frm; i++)
    {
        j_E_min[0] = - j_E[i][0], j_E_min[1] = - j_E[i][1], j_E_min[2] = - j_E[i][2];
        K_rot_rot(j_E_min, x, cells[i], norm[i], num_param, ed_param.k[1], cur_res3);
        sum_K[0] += cur_res3[0];
        sum_K[1] += cur_res3[1];
        sum_K[2] += cur_res3[2];


        j_M_min[0] = - j_M[i][0], j_M_min[1] = - j_M[i][1], j_M_min[2] = - j_M[i][2];
        R_rot_Far(j_M_min, x, cells[i], num_param, ed_param.k[1], cur_res3);
        sum_R[0] += cur_res3[0];
        sum_R[1] += cur_res3[1];
        sum_R[2] += cur_res3[2];
    }
    field_E[0] = deg * sum_K[0] - sum_R[0];
    field_E[1] = deg * sum_K[1] - sum_R[1];
    field_E[2] = deg * sum_K[2] - sum_R[2];
}

#endif