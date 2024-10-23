#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "common_type.h"
#include "get_j.h"
#include "e0.h"
#include "get_field.h"
#include "ED_Par.h"
#include "Num_Par.h"
#include "write_to_file.h"

int main(int argc, char **argv)
{
    const std::string input_dir = argc > 1 ? argv[1] : "../data/", \
                      result_dir = argc > 2 ? argv[2] : "../results/", \
                      geom_fname = argc > 3 ? argv[3] : "geodat_40_40.dat", \
                      ed_param_fname = argc > 4 ? argv[4] : "ed_param.txt", \
                      num_parFar_fname = argc > 5 ? argv[5] : "num_parFar_field.txt", \
                      num_parNear_fname = argc > 6 ? argv[6] : "num_parNear_field.txt";
    auto start = std::chrono::high_resolution_clock::now();

    //==========================Геометрия==============================
    TGrid_DC_Full a(input_dir + geom_fname);
    std::cout << "grid step: " << a.grid_step << std::endl;

    //=========================Parameters==============================
    ED_Par ed_param(input_dir + ed_param_fname);
    //E0
    E0_PolarType e0_H(0, ed_param.k_vec);
    E0_PolarType e0_V(1, ed_param.k_vec);

    //=======================Сохранение ed_Par=========================
    write_ED_Params_toFile(result_dir + "ed_param_res.txt", ed_param);

    //========================Задания токов H,V========================
    std::vector<std::complex<double>[3]> j_vec_H(a.num_frm), j_vec_V(a.num_frm);
    get_j_from_files(input_dir + "j/j_H_real.gv", input_dir + "j/j_H_image.gv", j_vec_H);
    get_j_from_files(input_dir + "j/j_V_real.gv", input_dir + "j/j_V_image.gv", j_vec_V);


    //===================Задания ячеек и нормалей =====================
    // Создадим список ячеек: cells[num_frm][4][3] и заполним его из сетки
    std::vector<double[4][3]> cells(a.num_frm);
    std::vector<double[3]> norm(a.num_frm);
    for (int i = 0; i < a.num_frm; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cells[i][j][0] = a.cell_list[i].cell[j][0];
            cells[i][j][1] = a.cell_list[i].cell[j][1];
            cells[i][j][2] = a.cell_list[i].cell[j][2];  
        }
        norm[i][0] = a.cell_list[i].norm[0];
        norm[i][1] = a.cell_list[i].norm[1];
        norm[i][2] = a.cell_list[i].norm[2];
    }

    //=====================Запись точек коллокации в файл===========
    std::string filename = result_dir + "body.gr";
    write_rkt_toFile(filename, cells);

    //=================Точки, в которых считаем поле===================
    //Открытие файла с точками(для примера те же, что и rkt)
    std::ifstream fin(result_dir + "body.gr");
    if (!fin.is_open()) {
        throw std::runtime_error("Can't opent file: " + result_dir + "body.gr");
    }

    int n1, n2, num_points;
    fin >> n1 >> n2;
    num_points = n1 * n2;//количество точек, в которых хотим считать поле 

    std::vector<double[3]> points_for_field(num_points);
    for (int i = 0; i < num_points; i++)
    {
        fin >> points_for_field[i][0] >> points_for_field[i][1] >> points_for_field[i][2];
    }

    fin.close();


    //===============Нахождение поля в заданных точках=================
    std::vector<std::complex<double>[3]> field_H(num_points), field_V(num_points);
    Num_Par num_parFar(input_dir + num_parFar_fname);
    Num_Par num_parNear(input_dir + num_parNear_fname);

    get_field_idealColloc(cells, norm, j_vec_H, ed_param.k[0], ed_param.k_vec, e0_H.e0,
                num_parFar, num_parNear, a.grid_step, points_for_field, field_H);
    get_field_idealColloc(cells, norm, j_vec_V, ed_param.k[0], ed_param.k_vec, e0_V.e0,
                num_parFar, num_parNear, a.grid_step, points_for_field, field_V);

    //====================Сохранение в формате .gv=====================
    filename = result_dir + "field/E_H";
    write_field_toFiles(filename, field_H, n1, n2);
    filename = result_dir + "field/E_V";
    write_field_toFiles(filename, field_V, n1, n2);


    //=========================Очистка памяти===================
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}

