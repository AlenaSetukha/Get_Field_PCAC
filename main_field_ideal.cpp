#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <chrono>
#include <vector>


#include "get_j.h"
#include "common_type.h"
#include "e0.h"

#include "get_field.h"
#include "ED_Par.h"
#include "Num_Par.h"

int main(int argc, char **argv)
{
    const std::string input_dir = argc > 1 ? argv[1] : "./data/", \
                      result_dir = argc > 2 ? argv[2] : "./results", \
                      geom_fname = argc > 3 ? argv[3] : "geodat_40_40.dat";
    auto start = std::chrono::high_resolution_clock::now();
//======================Геометрия==============================
    TGrid_DC_Full a(input_dir + geom_fname);
    a.fill_TGrid(input_dir + geom_fname);
    std::cout << "max diag: " << a.max_diag << std::endl;

//======================Parameters=============================
    ED_Par ed_param(input_dir + "ed_param.txt");
    //E0
    E0_Polar_Type e0_H(0, ed_param.k_vec);
    E0_Polar_Type e0_V(1, ed_param.k_vec);

//=====================Запись точек коллокации в файл===========
    std::ofstream fout(result_dir + "body.gr");
    if (!fout.is_open()) {
        throw std::runtime_error("Can't opent file: " + result_dir + "body.gr");
    }
    fout << 2 << " " << a.num_frm / 2 << std::endl;

    for (int i = 0; i < a.num_frm; i++)
    {
        fout << a.cell_list[i].rkt[0] << " " << a.cell_list[i].rkt[1] << " " << a.cell_list[i].rkt[2] << std::endl;
    }
    fout.close();

//========================Задания токов H,V========================
    std::vector<std::complex<double>[3]> j_vec_H(a.num_frm), j_vec_V(a.num_frm);
    get_j_from_files(input_dir + "j//j_H_real.gv", input_dir + "/j/j_H_image.gv", j_vec_H);
    get_j_from_files(input_dir + "j//j_V_real.gv", input_dir + "/j/j_V_image.gv", j_vec_V);


//===================Задания ячеек и нормалей =====================
    Num_Par num_param_field(input_dir + "/num_param_field.txt");
    
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

//=================Точки, в которых считаем поле===================
    //Открытие файла с точками
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

    get_field_ideal(cells, norm, j_vec_H, a.max_diag, a.num_frm,
                ed_param, num_param_field, e0_H.e0, points_for_field, field_H);
    get_field_ideal(cells, norm, j_vec_V, a.max_diag, a.num_frm,
                ed_param, num_param_field, e0_V.e0, points_for_field, field_V);

//====================Сохранение в формате .gv=====================
    save_field_as_gv(result_dir + "field/E_H", n1, n2, field_H);
    save_field_as_gv(result_dir + "field/E_V", n1, n2, field_V);

//=========================Очистка памяти===================


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}

