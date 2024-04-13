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

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
//======================Геометрия==============================
    TGrid_DC_Full a("data//geodat_40_40.dat");
    std::cout << "max diag: " << a.max_diag << std::endl;

//======================Parameters=============================
    ED_Par ed_param("data//ed_param.txt");
    //E0
    E0_Polar_Type e0_H(0, ed_param.k_vec);
    E0_Polar_Type e0_V(1, ed_param.k_vec);

//=====================Запись точек коллокации в файл===========
    std::ofstream fout("results//body.gr");
    if (!fout.is_open()) {
        std::cout << "Open results//body.gr error" << std::endl;
        return -1;
    }
    fout << 2 << " " << a.num_frm / 2 << std::endl;

    for (int i = 0; i < a.num_frm; i++)
    {
        fout << a.cell_list[i].rkt[0] << " " << a.cell_list[i].rkt[1] << " " << a.cell_list[i].rkt[2] << std::endl;
    }
    fout.close();

//========================Задания токов H,V========================
    std::vector<std::complex<double>[3]> j_vec_H(a.num_frm), j_vec_V(a.num_frm);
    std::string input_dir = "data//";

    get_j_from_files(input_dir + "j//j_H_real.gv", input_dir + "/j/j_H_image.gv", j_vec_H);
    get_j_from_files(input_dir + "j//j_V_real.gv", input_dir + "/j/j_V_image.gv", j_vec_V);

//====================Задания ячеек и нормалей=====================
    Num_Par num_param_field("data//num_param_field.txt");
    
    std::vector<double[4][3]> cells(a.num_frm);
    std::vector<double[3]> norm(a.num_frm);

    for (int i = 0; i < a.num_frm; i++)
    {
        norm[i][0] = a.cell_list[i].norm[0], norm[i][1] = a.cell_list[i].norm[1], norm[i][2] = a.cell_list[i].norm[2];
        for (int j = 0; j < 4; j++)
        {
            cells[i][j][0] = a.cell_list[i].cell[j][0];
            cells[i][j][1] = a.cell_list[i].cell[j][1];
            cells[i][j][2] = a.cell_list[i].cell[j][2];
        }
    }

//=================Точки, в которых считаем поле===================
    //Открытие файла с точками
    std::ifstream fin("results//body.gr");
    if (!fin.is_open()) {
        std::cout << "Open  results//body.gr  error" << std::endl;
        return -1;
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

    get_field_ideal(cells, norm, j_vec_H, a.max_diag, ed_param, num_param_field,
            e0_H.e0, points_for_field, field_H);
    get_field_ideal(cells, norm, j_vec_V, a.max_diag, ed_param, num_param_field,
            e0_V.e0, points_for_field, field_V);

//====================Сохранение в формате .gv=====================
    save_field_as_gv("results//field//E_H", n1, n2, field_H);
    save_field_as_gv("results//field//E_V", n1, n2, field_V);

//=========================Очистка памяти===================
    
    norm.clear();
    cells.clear();
    j_vec_H.clear();
    j_vec_V.clear();
    points_for_field.clear();
    field_H.clear();
    field_V.clear();




    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Программа выполнялась: " << duration.count() << " секунд" << std::endl;
    return 0;
}

