#ifndef _GET_FIELD_H_
#define _GET_FIELD_H_

#include <complex>
#include <string>
#include <vector>

#include "K.h"
#include "element_geom.h"
#include "integral_par.h"
#include "f_par.h"
#include "get_einc.h"
#include "ED_Par.h"
#include "Num_Par.h"


//===============================================================================================
//-----------------Функция расчета поля в точках из файла для идеального провдоника--------------
//===============================================================================================
//---------------------E(x) = sumj(K[sigma[j]], j_vec[j](x)) + Einc(x)---------------------------
// cells - координаты четырехугольных ячеек: [num_frm][4][3]
// norm - нормали к каждой ячейке: [num_frm][3]
// j_vec - вектор с электрическими токами: [num_frm][3]
// max_diag - шаг сетки(важен для сглаживания ядер)
// num_frm - число ячеек
// ed_param - э/д параметры
// num_param - численные параметры для инетгралов
// e0 - вектор поляризации 
// n_points - число точек для рассчета поля
// points_for_field - точки, в которых считается поле: [n_points][3]
// field - комплексное поле в каждой точке(результат): [n_points][3]

template<size_t CellPoints>
void get_field_ideal(const std::vector<double[CellPoints][3]> &cells, const std::vector<double[3]> &norm,
        const std::complex<double>** j_vec,
        const double max_diag, const int num_frm,
        const ED_Par& ed_param, const Num_Par& num_param,
        const double* e0,
        const std::vector<double[3]> &points_for_field,
        std::vector<std::complex<double>[3]> &field)
{
    size_t n_points = points_for_field.size();
    //Инициализация параметров
    //f_simple_pot_G
    integral_par f_simple_pot_G_par(1, num_param.n_start, num_param.p_max, num_param.eps);//1,10,1,0.001
    f_par param(max_diag * num_param.rs, ed_param.k[0]);

    //f_grad_simple_pot_G
    integral_par f_grad_simple_pot_G_par(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);//3,10,1,0.001
    f_par param_seg(max_diag * num_param.rs, ed_param.k[0]);

//========================Вычисление поля===========================

#pragma omp parallel for
    for (size_t i = 0; i < n_points; i++)
    {   
        std::complex<double> u[3], cur_res3[3];
        double cell_j[CellPoints][3], norm_j[3];
        std::complex<double> j_vec_j[3];
        field[i][0] = 0., field[i][1] = 0., field[i][2] = 0.;
        for (int j = 0; j < num_frm; j++)
        {
            for (int g = 0; g < CellPoints; g++)
            {
                cell_j[g][0] = cells[j][g][0];
                cell_j[g][1] = cells[j][g][1];
                cell_j[g][2] = cells[j][g][2];
            }
            norm_j[0] = norm[j][0], norm_j[1] = norm[j][1], norm_j[2] = norm[j][2];
            j_vec_j[0] = j_vec[j][0], j_vec_j[1] = j_vec[j][1], j_vec_j[2] = j_vec[j][2];
            K_rot_rot(j_vec_j, points_for_field[i], cell_j, norm_j, f_simple_pot_G_par,
                    f_grad_simple_pot_G_par, param, param_seg, cur_res3);
            field[i][0] += cur_res3[0];
            field[i][1] += cur_res3[1];
            field[i][2] += cur_res3[2];
        }
        //вычисление e_inc
        get_einc(points_for_field[i], e0, ed_param.k_vec, u);
        field[i][0] += u[0];
        field[i][1] += u[1];
        field[i][2] += u[2];
    }
    return;
}




void save_field_as_gv(const std::string filename_out,
    const int n1, const int n2,
    std::vector<std::complex<double>[3]> &field)
{
    //Создание файлов для записи электрического поля
    std::ofstream fout_u_real(filename_out + "_real.gv");              //Электрическое поле, действительная часть
    if (!fout_u_real.is_open()) {
        std::cout << "Open " + filename_out + "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_im(filename_out + "_im.gv");                  //Электрическое поле, мнимая часть
    if (!fout_u_im.is_open()) {
        std::cout << "Open " + filename_out + "_im.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_abs(filename_out + "_abs.gdr");               //Модуль электрического поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " + filename_out + "_abs.gdr error" << std::endl;
        return;
    }

    // Запись
    double u_abs;
    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_im << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    int n = n1 * n2;
    for (int i = 0; i < n; i++)
    {
        fout_u_real << field[i][0].real() << " " << field[i][1].real() << " " << field[i][2].real() << std::endl;
        fout_u_im << field[i][0].imag() << " " << field[i][1].imag() << " " << field[i][2].imag() << std::endl;
        u_abs = sqrt(abs_tmp(field[i][0]) * abs_tmp(field[i][0]) +
                    abs_tmp(field[i][1]) * abs_tmp(field[i][1])  +
                    abs_tmp(field[i][2]) * abs_tmp(field[i][2]));
        fout_u_abs << u_abs << std::endl;
    }

    // Закрытие файлов
    fout_u_real.close();
    fout_u_im.close();
    fout_u_abs.close();
    return;
}








/*
//===============================================================================================
//-----------------Функция расчета поля в точках из файла для идеального провдоника--------------
//===============================================================================================
// a - класс, описывающий сетку
// filename_in - имя файла с точками расчета
// j_vec - вектор с электрическими токами
// ed_param - э/д параметры(каждой области)
// num_param - численные параметры для инетгралов
// e0 - вектор поляризации 
// filename_out - файл, в который записываем результат

template<typename TGrid>
void get_field_ideal(const TGrid& a, const std::string filename_in, 
        const std::complex<double>** j_vec,
        const ED_Par& ed_param, const Num_Par& num_param,
        const double* e0,
        const std::string filename_out)
{
    //Инициализация параметров
    //f_simple_pot_G
    integral_par f_simple_pot_G_par(1, num_param.n_start, num_param.p_max, num_param.eps);//1,10,1,0.001

    f_par param(a.max_diag * num_param.rs, ed_param.k[0]);

    //f_grad_simple_pot_G
    integral_par f_grad_simple_pot_G_par(3, num_param.n_start_seg, num_param.p_max_seg, num_param.eps);//3,10,1,0.001

    //f_par param_seg(0.000000000001, ed_param.k[0]);
    f_par param_seg(a.max_diag * num_param.rs, ed_param.k[0]);



    //Создание файлов для записи электрического поля
    std::ofstream fout_u_real(filename_out + "_real.gv");              //Электрическое поле, действительная часть
    if (!fout_u_real.is_open()) {
        std::cout << "Open " + filename_out + "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_im(filename_out + "_im.gv");                  //Электрическое поле, мнимая часть
    if (!fout_u_im.is_open()) {
        std::cout << "Open " + filename_out + "_im.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_abs(filename_out + "_abs.gdr");               //Модуль электрического поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " + filename_out + "_abs.gdr error" << std::endl;
        return;
    }

    //Открытие файла с точками
    std::ifstream fin(filename_in);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " + filename_in + " error" << std::endl;
        return;
    }



//========================Вычисление поля===========================
    Constants c;
    std::complex<double> deg, deg1;

    int n1, n2, n;
    fin >> n1 >> n2;

    n = n1 * n2;//количество точек, в которых хотим считать поле 
    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_im << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;



    double x_0[3]; // точка, в которой считаем поле
    std::complex<double> u[3];
    std::complex<double> cur_res3[3], j_vec_j[3];
    double u_abs;

    for (int i = 0; i < n; i++)
    {
        fin >> x_0[0] >> x_0[1] >> x_0[2];
        //вычисление e_inc, h_inc
        get_einc(x_0, e0, ed_param.k_vec, u);

        for (int j = 0; j < a.num_frm; j++)
        {
            j_vec_j[0] = j_vec[j][0], j_vec_j[1] = j_vec[j][1], j_vec_j[2] = j_vec[j][2];
            K_rot_rot(j_vec_j, x_0, a.cell_list[j].cell, a.cell_list[j].norm,
                        f_simple_pot_G_par, f_grad_simple_pot_G_par, param, param_seg, cur_res3);
            u[0] += cur_res3[0];
            u[1] += cur_res3[1];
            u[2] += cur_res3[2];
        }

        fout_u_real << u[0].real() << " " << u[1].real() << " " << u[2].real() << std::endl;
        fout_u_im << u[0].imag() << " " << u[1].imag() << " " << u[2].imag() << std::endl;
        u_abs = sqrt(abs_tmp(u[0]) * abs_tmp(u[0]) + abs_tmp(u[1]) * abs_tmp(u[1])  + abs_tmp(u[2]) * abs_tmp(u[2]));
        fout_u_abs << u_abs << std::endl;
    }

    fin.close();
    fout_u_real.close();
    fout_u_im.close();
    fout_u_abs.close();
    return;
}





#include "get_area.h"
#include "get_field_one_point.h"
//===============================================================================================
//----------------------------Функции расчета электрического поля для диэлектрика----------------
//===============================================================================================
// a - класс, описывающий сетку
// filename_in - имя файла с точками расчета
// j_E - вектор с электрическими токами
// j_M - вектор с магнитными токами
// ed_param - э/д параметры(каждой области)
// num_param - численные параметры для инетгралов
// e0 - вектор поляризации 
// filename_out - файл, в который записываем результат

template<typename TGrid>
void get_field_dielectric(const TGrid& a, const std::string filename_in,
        const std::complex<double>** j_E , const std::complex<double>** j_M,
        const ED_Par& ed_param, const Num_Par& num_param,
        const double* e0,
        const std::string filename_out)
{
    //Создание файлов с результатами
    std::ofstream fout_u_real(filename_out + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_out << "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_image(filename_out + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_out << "_image.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_abs(filename_out + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_out << "_abs.gdr error" << std::endl;
        return;
    }



    std::ifstream fin(filename_in);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename_in << " error" << std::endl;
        return;
    }

    //Расчет поля в точках из файла
    int n, n1, n2, area;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3];

    for (int i = 0; i < n; i++)
    {
        fin >> x[0] >> x[1] >> x[2];
        //Определение, в какой области лежит точка(вне или внутри диэлектрика)
        area = get_area(x, a);
        if (area == -1)
        {
            get_field_dielectric_area2(a, x, j_E, j_M, ed_param, num_param, field_val);
        } else {
            get_field_dielectric_area1(a, x, j_E, j_M, ed_param, num_param, e0, field_val);
        }
        fout_u_real << field_val[0].real() << " " << field_val[1].real() << " " << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() << " " << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }

    //Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return;
}









//===============================================================================================
//----------------------------Дополнительные функции для диэлектрика-----------------------------
//===============================================================================================

//------------------------------Поле внутри = Einc, снаружи = 0----------------------------------
template<typename TGrid>
void get_field_exact(const TGrid& a, const std::string filename_in,
        const ED_Par& ed_param, const double* e0,
        const std::string filename_out)
{
    // Создание файлов
    std::ofstream fout_u_real(filename_out + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_out << "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_image(filename_out + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_out << "_image.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_abs(filename_out + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_out << "_abs.gdr error" << std::endl;
        return;
    }

    std::ifstream fin(filename_in);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename_in << " error" << std::endl;
        return;
    }


    // Расчет поля в точках из файла
    int n, n1, n2, area;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3], deg;

    for (int i = 0; i < n; i++)
    {
        fin >> x[0] >> x[1] >> x[2];
        //Определение, в какой области лежит точка
        area = get_area(x, a);
        if (area == 1)//поле снаружи = 0
        {
            field_val[0] = std::complex<double>(0., 0.);
            field_val[1] = std::complex<double>(0., 0.);
            field_val[2] = std::complex<double>(0., 0.);
        } else
        {
            get_einc(x, e0, ed_param.k_vec, field_val);
        }
        fout_u_real << field_val[0].real() << " " << field_val[1].real() <<
                        " " << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() <<
                        " " << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }

    // Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return;
}








//---------------Поле в любой точке диэлектрика равно полю во внутренней области-----------------
template<typename TGrid>
void get_field_omega2(const TGrid& a, const std::string filename_in,
        const std::complex<double>** j_E , const std::complex<double>** j_M,
        const ED_Par& ed_param, const Num_Par& num_param,
        const std::string filename_out)
{
//=============================Создание файлов===================================
    std::ofstream fout_u_real(filename_out + "_real.gv");
    if (!fout_u_real.is_open()) {
        std::cout << "Open " << filename_out << "_real.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_image(filename_out + "_image.gv");
    if (!fout_u_image.is_open()) {
        std::cout << "Open " << filename_out << "_image.gv error" << std::endl;
        return;
    }

    std::ofstream fout_u_abs(filename_out + "_abs.gdr"); //Модуль электрическое поля
    if (!fout_u_abs.is_open()) {
        std::cout << "Open " << filename_out << "_abs.gdr error" << std::endl;
        return;
    }

    std::ifstream fin(filename_in);                // Точки, в которых считаем поле
    if (!fin.is_open()) {
        std::cout << "Open " << filename_in << " error" << std::endl;
        return;
    }

//========================Расчет поля в точках из файла===========================

    int n, n1, n2;
    fin >> n1 >> n2;


    fout_u_real << n1 << " " << n2 << std::endl;
    fout_u_image << n1 << " " << n2 << std::endl;
    fout_u_abs << n1 << " " << n2 << std::endl;

    Constants c;
    n = n1 * n2;//количество точек, в которых хотим считать поле
    double x[3]; // точка, в которой считаем поле
    std::complex<double> field_val[3];

    for (int i = 0; i < n; i++)
    {
        fin >> x[0] >> x[1] >> x[2];
        get_field_dielectric_area2(a, x, j_E, j_M, ed_param, num_param, field_val);
        fout_u_real << field_val[0].real() << " " << field_val[1].real() << " "
                        << field_val[2].real() << std::endl;
        fout_u_image << field_val[0].imag() << " " << field_val[1].imag() << " "
                        << field_val[2].imag() << std::endl;
        fout_u_abs << vec_length(field_val) << std::endl;
    }

    //Закрытие файлов
    fout_u_real.close();
    fout_u_image.close();
    fout_u_abs.close();
    fin.close();
    return;
}
*/
#endif
