#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#include "element_geom.h"
#include "common_type.h"
#include "cell.h"

// sum - число ячеек
// ssum - общее число отрезков на всех линиях отрыва
// kl - число линий отрыва
// kseg - число отрезков на одной из линий отрыва(tmp)
// t - счетчик для подсчета модулей на объекте

//===============================================================================================
//------------------------------------Constructor------------------------------------------------
//===============================================================================================
TGrid_DC_Full::TGrid_DC_Full(const std::string filename)
{
    max_diag = 0.;
    num_mod = 0;
    num_obj = 0;
    num_bl = 0;
    num_seg = 0;

    int t = 0, sum = 0;
    int kmod = 0;

    int sign = 0, nfm = 0, i11 = 0, i22 = 0;
    int kseg = 0;

    double root[4][3];
    Cell cell_tmp;

    std::ifstream fin2(filename);
    if (!fin2.is_open()) {
        std::cout << "Read geodat.dat error" << std::endl;
        exit(1);
    }


    // Заполнение сетки
    std::array<int, 2> be_i, beo_i;

    fin2 >> num_obj;
    for (int i = 0; i < num_obj; i++)
    {
        fin2 >> kmod;
        num_mod += kmod;
        beo_i[0] = t;
        for (int j = 0; j < kmod; j++)
        {
            fin2 >> sign >> nfm >> i11 >> i22;
            be_i[0] = sum;
            for (int k = 0; k < nfm; k++)
            {
                for (int g = 0; g < 4; g++)
                {
                    fin2 >> root[g][0] >> root[g][1] >> root[g][2];
                }
                cell_tmp.Cell_fill(root);
                cell_list.push_back(cell_tmp);
                //Нахождение шага сетки
                check_step(cell_list[sum].cell, this->max_diag);
                sum++;                                                             // number of frames
            }
            be_i[1] = sum - 1;
            beg_end.push_back(be_i);
            t++;
        }
        beo_i[1] = t - 1;
        beg_end_obj.push_back(beo_i);
    }
    num_frm = sum;

    // Заполнений линий отрыва
    std::array<std::array<double, 3>, 2> xb_i;
    fin2 >> num_bl;                                                                    // число линий отрыва
    for (int i = 0; i < num_bl; i++)
    {
        fin2 >> kseg;                                                                  // кол-во отрезков на каждой линии отрыва
        num_seg += kseg;                                                             
        for (int j = 0; j < kseg; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                fin2 >> xb_i[k][0] >> xb_i[k][1] >> xb_i[k][2];
            }
            x_bound.push_back(xb_i);
        }
    }
    fin2.close();
}




//===============================================================================================
//-----------------------------------Copy Constructor--------------------------------------------
//===============================================================================================
TGrid_DC_Full::TGrid_DC_Full(const TGrid_DC_Full& obj)
{
    num_obj = obj.num_obj;
    num_mod = obj.num_mod;
    num_frm = obj.num_frm;
    max_diag = obj.max_diag;

    num_bl = obj.num_bl;
    num_seg = obj.num_seg;

    for (int i = 0; i < num_frm; i++)
    {
        cell_list.push_back(obj.cell_list[i]);
    }
   
    for (int i = 0; i < num_seg; i++)
    {
        x_bound.push_back(obj.x_bound[i]);
    }

    for (int i = 0; i < num_mod; i++)
    {
        beg_end.push_back(obj.beg_end[i]);
    }

    for (int i = 0; i < num_obj; i++)
    {
        beg_end_obj.push_back(obj.beg_end_obj[i]);
    }
}




//===============================================================================================
//-----------------------------------Destructor--------------------------------------------------
//===============================================================================================
TGrid_DC_Full::~TGrid_DC_Full()
{
    cell_list.clear();
    x_bound.clear();
    beg_end.clear();
    beg_end_obj.clear();
}





//===============================================================================================
//-----------------------------------Step searching----------------------------------------------
//===============================================================================================
void TGrid_DC_Full::check_step(const double (&root_tmp)[4][3], double max_diag_tmp)
{
    double diam_cur = get_diam(root_tmp);
    if (diam_cur > max_diag_tmp)
    {
        this->max_diag = diam_cur;      //заполнение шага сетки
    }
}
