#ifndef _GET_FIELD_H_
#define _GET_FIELD_H_

#include <complex>
#include <vector>

#include "integral_par.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "get_fieldOnePoint.h"
#include "get_area.h"


//===================================================================================
//-----------Электрическое поле для идеального провдоника в дальней зоне-------------
//===================================================================================
/**
 * Расчет полного эл. поля только по электрическим токам в дальней зоне
 *                 E(x) = sumj(K[sigma[j]], j_vec[j](x)) + Einc(x)
 * K = K_far
 * Параметры функции:
 *      cells[num_frm][4/3][3] - координаты ячеек
 *      j_vec[num_frm][3] - векторы электрических токов в точках коллокации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      e0[3] - вектор внешней поляризации
 *      num_param - численные параметры для расчета интегралов
 *      x[n_pnt][3] - список точек для рассчета в них эл. поля
 *      field_E[n_pnt][3] - результат, вектор эл. поля в точках
 * 
 * Примечания:
 *      num_param.rs - радиус сглаживания (интеграл по ячейке) относительно диам. ячейки
 */
template<size_t CellPoints>
void get_field_idealFar(const std::vector<double[CellPoints][3]> &cells,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::complex<double> k0,
        const double* k_vec, const double* e0,
        const Num_Par& num_param,
        const std::vector<double[3]> &x,
        std::vector<std::complex<double>[3]> &field_E)
{
    //===========Вычисление поля========
    int n_points = x.size();
    for (int i = 0; i < n_points; i++)
    {
        get_fieldOnePoint_idealFar(x[i], e0, k0, k_vec, j_vec,
                cells, num_param, field_E[i]);
    }
}





//===================================================================================
//-------Электрическое поле для идеального провдоника в дальней зоне с уточнением----
//===================================================================================
/**
 * Полное эл. поле по электрическим токам в дальней зоне с уточнением
 *                 E(x) = sumj(K[sigma[j]], j_vec[j](x)) + Einc(x)
 * K = K_far или K
 * Параметры функции:
 *      cells[num_frm][4/3][3] - координаты ячеек
 *      j_vec[num_frm][3] - векторы электрических токов в точках коллокации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      e0[3] - вектор внешней поляризации
 *      num_parFar - численные параметры(для интегрирования) дальней зоны
 *      num_parNear - параметры ближней зоны без выделения особенности
 *      field_E[num_frm][3] - результат, поле в точке 
 * 
 * num_parFar.rs - радиус сглаживания (интеграл по ячейке)
 * относительно одной ячейки(~0.5 - 1)
 * 
 * num_parNear.rs, rs_seg - радиус сглаживания (интеграл по ячейке/отрезку)
 * относительно (diamj / n_start) (~0.5-2)
 */


template<size_t CellPoints>
void get_field_idealFarClar(const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::complex<double> k0,
        const double* k_vec, const double* e0,
        const Num_Par& num_parFar, const Num_Par& num_parNear,
        const double grid_step,
        const std::vector<double[3]> &x,
        std::vector<std::complex<double>[3]> &field_E)
{
    //===========Вычисление поля========
    int n_points = x.size();
    for (int i = 0; i < n_points; i++)
    {
        get_fieldOnePoint_idealFarClar(x[i], e0, k0, k_vec, j_vec,
                cells, norm, num_parFar, num_parNear, grid_step, field_E[i]);
    }
}






//===================================================================================
//---------Электрическое поле для идеального провдоника в точках коллокаций----------
//===================================================================================
/**
 * Полное эл. поле по электрическим токам в дальней зоне с уточнением
 *                 E(x) = sumj(K[sigma[j]], j_vec[j](x)) + Einc(x)
 * K = K_far или K_near(с выделением особенности)
 * Параметры функции:
 *      cells[num_frm][4/3][3] - координаты ячеек
 *      j_vec[num_frm][3] - векторы электрических токов в точках коллокации
 *      k0 - внешнее волновое число
 *      k_vec[3] - внешний волновой вектор
 *      e0[3] - вектор внешней поляризации
 *      num_parFar - численные параметры(для интегрирования) дальней зоны
 *      num_parNear - параметры ближней зоны без выделения особенности
 *      field_E[num_frm][3] - результат, поле в точке 
 * 
 * num_parFar.rs - радиус сглаживания (интеграл по ячейке)
 * относительно одной ячейки(~0.5 - 1)
 * 
 * num_parNear.rs, rs_seg - радиус сглаживания (интеграл по ячейке/отрезку)
 * относительно (diamj / n_start) (~0.5-2)
 */


template<size_t CellPoints>
void get_field_idealColloc(const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const std::vector<std::complex<double>[3]> &j_vec,
        const std::complex<double> k0,
        const double* k_vec, const double* e0,
        const Num_Par& num_parFar, const Num_Par& num_parNear,
        const double grid_step,
        const std::vector<double[3]> &x,
        std::vector<std::complex<double>[3]> &field_E)
{
    //===========Вычисление поля========
    int n_points = x.size();
    for (int i = 0; i < n_points; i++)
    {
        get_fieldOnePoint_idealColloc(x[i], e0, k0, k_vec, j_vec,
                cells, norm, num_parFar, num_parNear, grid_step, field_E[i]);
    }
}










//===================================================================================
//------------Электрическое поле в точках вокруг одного диэлектрика------------------
//===================================================================================
/**
 * Расчет полного эл. поля для одного диэлектрика
 * Параметры функции:
 *      cells[num_frm][4/3][3] - координаты ячеек
 *      norm[num_frm][3] - векторы нормалей к ячейкам
 *      j_E/j_M[num_frm][3] - векторы э/м токов в точках коллокации
 *      e0[3] - вектор внешней поляризации
 *      ed_param - э/д параметры задачи
 *      num_param - ичсленные параметры(для интегрирования)
 *      field_E[num_frm][3] - результат, поле в точке 
 * 
 * 
 * num_param.rs - сглаживание относительно каждой ячейки(~0.5-1)
 */

template<size_t CellPoints>
void get_fieldDielectric(const std::vector<double[CellPoints][3]> &cells,
        const std::vector<double[3]> &norm,
        const std::vector<std::complex<double>[3]> &j_E,
        const std::vector<std::complex<double>[3]> &j_M,
        const double* e0,
        const ED_Par& ed_param, const Num_Par& num_param,
        const std::vector<double[3]> &x,
        std::vector<std::complex<double>[3]> &field_E)
{
    int area, n_points = x.size();
    for (int i = 0; i < n_points; i++)
    {
        // Вне или внутри диэлектрика
        area = get_area(x[i], cells, norm);
        if (area == -1)
        {
            get_fieldDielectricA1(x[i], e0, cells, norm,
                    j_E, j_M, ed_param, num_param, field_E[i]);
        } else {
            get_fieldDielectricA2(x[i], cells, norm,
                    j_E, j_M, ed_param, num_param, field_E[i]);
        }
    }
}
#endif