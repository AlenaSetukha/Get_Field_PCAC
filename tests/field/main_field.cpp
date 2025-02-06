#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>
#include <omp.h>
#include <cblas.h>
#include <lapacke.h>
#include <openblas_config.h>



#include "common_type.h"
#include "E0.h"
#include "get_field.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "write_to_file.h"
#include "get_einc.h"
#include "get_hinc.h"


//Тестовый расчет электромагнитного поля 


int main(int argc, char **argv)
{
    const std::string input_dir = argc > 1 ? argv[1] : "../tests/data/", \
                      result_dir = argc > 2 ? argv[2] : "../tests/results/", \
                      geom_fname = argc > 3 ? argv[3] : "grids/30_50_1_1_1.dat", \
                      ed_param_fname = argc > 4 ? argv[4] : "ed_param.txt";
    auto start = std::chrono::high_resolution_clock::now();

    //==========================Геометрия==============================
    TGrid_DC_Full a(input_dir + geom_fname);
    std::cout << "grid step: " << a.grid_step << std::endl;


    //=====================Запись точек коллокации в файл===========
    std::string filename = result_dir + "body.gr";
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Open " + filename + " error: write_rkt_toFile");
    }

    fout << 2 << " " << a.num_frm / 2 << std::endl;

    for (int i = 0; i < a.num_frm; i++) {
        fout << a.cell_list[i].rkt[0] << " " <<
                a.cell_list[i].rkt[1] << " " <<
                a.cell_list[i].rkt[2] << std::endl;
    }
    fout.close();

    //=========================Parameters==============================
    ED_Par ed_param(input_dir + ed_param_fname);
    E0_PolarType e0_H(0, ed_param.k_vec);

    //=======================Сохранение ed_Par=========================
    write_ED_Params_toFile(result_dir + "ed_param_res.txt", ed_param);




    //=============Задания токов на ячейках разбиения==================
    std::vector<std::complex<double>[3]> j_E(a.num_frm), j_M(a.num_frm);
    std::complex<double> H_inc[3], E_inc[3], vp[3];

    for (int i = 0; i < a.num_frm; i++) {
        get_hinc(a.cell_list[i].rkt, e0_H.e0, ed_param, H_inc);
        vec_prod(a.cell_list[i].norm, H_inc, vp);
        j_E[i][0] = -vp[0], j_E[i][1] = -vp[1], j_E[i][2] = -vp[2];


        get_einc(a.cell_list[i].rkt, e0_H.e0, ed_param.k_vec, E_inc);
        vec_prod(a.cell_list[i].norm, E_inc, vp);
        j_M[i][0] = vp[0], j_M[i][1] = vp[1], j_M[i][2] = vp[2];
    }
    







    //=================Точки, в которых считаем поле===================
    std::ifstream fin(input_dir + "grids/grid_OXY_4_4.gr");
    if (!fin.is_open()) {
        throw std::runtime_error("Can't opent file: " + input_dir + "grid.gr: main");
    }

    int n1, n2, num_points;
    fin >> n1 >> n2;
    num_points = n1 * n2;//количество точек, в которых хотим считать поле 

    std::vector<double[3]> points_for_field(num_points);
    for (int i = 0; i < num_points; i++) {
        fin >> points_for_field[i][0] >> points_for_field[i][1] >> points_for_field[i][2];
    }

    fin.close();


    //=====================Нахождение поля в лоб=======================
    std::vector<std::complex<double>[3]> field_E(num_points), field_H(num_points);
    for (int i = 0; i < num_points; i++) {
        get_einc(points_for_field[i], e0_H.e0, ed_param.k_vec, field_E[i]);
        get_hinc(points_for_field[i], e0_H.e0, ed_param, field_H[i]);
    }


    //==================Нахождение поля, Стреттон-Чу===================
    std::vector<std::complex<double>[3]> field_E_SCH(num_points), field_H_SCH(num_points);


    std::complex<double> sum_KE[3]{}, sum_RE[3]{};
    std::complex<double> sum_KH[3]{}, sum_RH[3]{};
    std::complex<double> cur_res3[3]{};
    std::complex<double> degE = std::complex<double>(0., 1.) /
                    (ed_param.omega0 * ed_param.eps_d[0] * Constants::eps0);

    std::complex<double> degH = std::complex<double>(0., 1.) /
                    (ed_param.omega0 * ed_param.mu_d[0] * Constants::m0);


    Num_Par num_param(0.00001, 1.5, 1.5, 8, 10, 2, 2);

    for (int i = 0; i < num_points; i++) {
        sum_KE[0] = 0., sum_KE[1] = 0., sum_KE[2] = 0.;
        sum_RE[0] = 0., sum_RE[1] = 0., sum_RE[2] = 0.;
        sum_KH[0] = 0., sum_KH[1] = 0., sum_KH[2] = 0.;
        sum_RH[0] = 0., sum_RH[1] = 0., sum_RH[2] = 0.;

        for (int j = 0; j < a.num_frm; j++) {
            // Поле E
            double pnt_dist = dist(points_for_field[i], a.cell_list[j].rkt);
            if (pnt_dist > 3. * a.grid_step) { // Если точка далеко от поверхности
                K_rot_rot_Far(j_E[j], points_for_field[i], a.cell_list[j].cell,
                                num_param, ed_param.k[0], cur_res3);
                sum_KE[0] += cur_res3[0];
                sum_KE[1] += cur_res3[1];
                sum_KE[2] += cur_res3[2];

                R_rot(j_M[j], points_for_field[i], a.cell_list[j].cell, num_param, ed_param.k[0], cur_res3);
                sum_RE[0] += cur_res3[0];
                sum_RE[1] += cur_res3[1];
                sum_RE[2] += cur_res3[2];

                // Поле H
                K_rot_rot_Far(j_M[j], points_for_field[i], a.cell_list[j].cell,
                                num_param, ed_param.k[0], cur_res3);
                sum_KH[0] += cur_res3[0];
                sum_KH[1] += cur_res3[1];
                sum_KH[2] += cur_res3[2];

                R_rot(j_E[j], points_for_field[i], a.cell_list[j].cell, num_param, ed_param.k[0], cur_res3);
                sum_RH[0] += cur_res3[0];
                sum_RH[1] += cur_res3[1];
                sum_RH[2] += cur_res3[2];
            } else { // если точка близко к поверхности
                K_rot_rot_Near(j_E[j], points_for_field[i], a.cell_list[j].cell,
                                num_param, ed_param.k[0], cur_res3);
                sum_KE[0] += cur_res3[0];
                sum_KE[1] += cur_res3[1];
                sum_KE[2] += cur_res3[2];

                R_rot(j_M[j], points_for_field[i], a.cell_list[j].cell, num_param, ed_param.k[0], cur_res3);
                sum_RE[0] += cur_res3[0];
                sum_RE[1] += cur_res3[1];
                sum_RE[2] += cur_res3[2];

                // Поле H
                K_rot_rot_Near(j_M[j], points_for_field[i], a.cell_list[j].cell,
                                num_param, ed_param.k[0], cur_res3);
                sum_KH[0] += cur_res3[0];
                sum_KH[1] += cur_res3[1];
                sum_KH[2] += cur_res3[2];

                R_rot(j_E[j], points_for_field[i], a.cell_list[j].cell, num_param, ed_param.k[0], cur_res3);
                sum_RH[0] += cur_res3[0];
                sum_RH[1] += cur_res3[1];
                sum_RH[2] += cur_res3[2];
            }
        }
        field_E_SCH[i][0] = degE * sum_KE[0] - sum_RE[0];
        field_E_SCH[i][1] = degE * sum_KE[1] - sum_RE[1];
        field_E_SCH[i][2] = degE * sum_KE[2] - sum_RE[2];


        field_H_SCH[i][0] = degH * sum_KH[0] + sum_RH[0];
        field_H_SCH[i][1] = degH * sum_KH[1] + sum_RH[1];
        field_H_SCH[i][2] = degH * sum_KH[2] + sum_RH[2];
    }





    //====================Сохранение в формате .gv=====================
    filename = result_dir + "field/E";
    write_field_toFiles(filename, field_E, n1, n2);
    filename = result_dir + "field/H";
    write_field_toFiles(filename, field_H, n1, n2);

    filename = result_dir + "field/E_SCH";
    write_field_toFiles(filename, field_E_SCH, n1, n2);
    filename = result_dir + "field/H_SCH";
    write_field_toFiles(filename, field_H_SCH, n1, n2);


    //=========================Очистка памяти===================
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}