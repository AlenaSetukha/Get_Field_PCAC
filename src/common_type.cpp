#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "element_geom.h"
#include "common_type.h"
#include "constants.h"
#include "cell.h"

// sum - число ячеек
// ssum - общее число отрезков на всех линиях отрыва
// kl - число линий отрыва
// kseg - число отрезков на одной из линий отрыва(tmp)
// t - счетчик для подсчета модулей на объекте



//===============================================================================================
//------------------------------------Constructor------------------------------------------------
//===============================================================================================
TGrid_DC_Full::TGrid_DC_Full(const std::string filename)      //сквозная сетка объектов и модулей
{
    int smod = 0, sum = 0, kobj = 0, kmod = 0, kl = 0, kseg = 0, sseg = 0;
    int  sign, nfm, i11, i22;
    double x;

//------------------------Calculation of the number of frames--------------------
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cout << "Read geodat.dat error" << std::endl;
        exit(-1);
    }

    fin >> kobj;
    for (int i = 0; i < kobj; i++)
    {
        fin >> kmod;
        smod += kmod;
        for (int j = 0; j < kmod; j++)
        {
            fin >> sign >> nfm >> i11 >> i22;                               // not used
            sum = sum + nfm;                                                // число ячеек
            nfm = nfm * 12;
            for (int k = 0; k < nfm; k++)
            {
                fin >> x;
            }
        }
    }

    fin >> kl;                                                              //число линий отрыва
    for (int i = 0; i < kobj; i++)
    {
        fin >> kseg;                                                        //кол-во отрезков на каждой линии отрыва
        sseg += kseg;                                                       //всего отрезков
        for (int j = 0; j < kseg; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                fin >> x;                                                   //not used
            }
        }
    }

    fin.close();
//-------------------------------Memory allocation-------------------------------
    num_obj = kobj;
    num_mod = smod;
    num_frm = sum;
    num_bl = kl;
    num_seg = sseg;


    cell_list = new Cell[num_frm];

    x_bound = new double**[num_seg];
    for (int i = 0; i < num_seg; i++) {
        x_bound[i] = new double*[2];
        for (int j = 0; j < 2; j++) {
            x_bound[i][j] = new double[3];
        }
    }

    beg_end = new int*[num_mod];
    for (int i = 0; i < num_mod; i++) {
        beg_end[i] = new int[2];
    }

    beg_end_obj = new int*[num_obj];
    for (int i = 0; i < num_obj; i++) {
        beg_end_obj[i] = new int[2];
    }

    return;
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

    cell_list = new Cell[num_frm];

    for (int i = 0; i < num_frm; i++)
    {
        cell_list[i].s = obj.cell_list[i].s;
        for (int j = 0; j < 3; j++)
        {
            cell_list[i].norm[j] = obj.cell_list[i].norm[j]; 
            cell_list[i].rkt[j] = obj.cell_list[i].rkt[j];
            cell_list[i].tau[0][j] = obj.cell_list[i].tau[0][j];
            cell_list[i].tau[1][j] = obj.cell_list[i].tau[1][j];
        }

        for (int j = 0; j < 4; j++)
        {
            cell_list[i].cell[j][0] = obj.cell_list[i].cell[j][0];
            cell_list[i].cell[j][1] = obj.cell_list[i].cell[j][1];
            cell_list[i].cell[j][2] = obj.cell_list[i].cell[j][2];
        }
    }
   
    x_bound = new double**[num_seg];
    for (int i = 0; i < num_seg; i++) {
        x_bound[i] = new double*[2];
        for (int j = 0; j < 2; j++) {
            x_bound[i][j] = new double[3];
            for (int k = 0; k < 3; k++)
            {
                x_bound[i][j][k] = obj.x_bound[i][j][k];
            }
        }
    }

    beg_end = new int*[num_mod];
    for (int i = 0; i < num_mod; i++) {
        beg_end[i] = new int[2];
        beg_end[i][0] = obj.beg_end[i][0];
        beg_end[i][1] = obj.beg_end[i][1];
    }

    beg_end_obj = new int*[num_obj];
    for (int i = 0; i < num_obj; i++) {
        beg_end_obj[i] = new int[2];
        beg_end_obj[i][0] = obj.beg_end_obj[i][0];
        beg_end_obj[i][1] = obj.beg_end_obj[i][1];
    }

    return;
}




//===============================================================================================
//-----------------------------------Destructor--------------------------------------------------
//===============================================================================================
TGrid_DC_Full::~TGrid_DC_Full()
{
    delete[] cell_list;

    for (int i = 0; i < num_seg; i++) {
        for (int j = 0; j < 2; j++) {
            delete[] x_bound[i][j];
        }
        delete[] x_bound[i];
    }
    delete[] x_bound;

    for (int i = 0; i < num_mod; i++) {
        delete[] beg_end[i];
    }
    delete[] beg_end;

    for (int i = 0; i < num_obj; i++) {
        delete[] beg_end_obj[i];
    }
    delete[] beg_end_obj;
}






//===============================================================================================
//---------------------------------Filling TGrid-------------------------------------------------
//===============================================================================================
void TGrid_DC_Full::fill_TGrid(const std::string filename)
{
    max_diag = 0.;
    int t = 0, sum = 0;
    int kobj = 0, kmod = 0, sign = 0, nfm = 0, i11 = 0, i22 = 0;
    int kl = 0, kseg = 0, ssum = 0;
    double root[4][3];

    std::ifstream fin2(filename);
    if (!fin2.is_open()) {
        std::cout << "Read geodat.dat error" << std::endl;
        exit(1);
    }

    fin2 >> kobj;
    for (int i = 0; i < kobj; i++)
    {
        fin2 >> kmod;
        beg_end_obj[i][0] = t;
        for (int j = 0; j < kmod; j++)
        {
            fin2 >> sign >> nfm >> i11 >> i22;
            beg_end[t][0] = sum;
            for (int k = 0; k < nfm; k++)
            {
                for (int g = 0; g < 4; g++)
                {
                    fin2 >> root[g][0] >> root[g][1] >> root[g][2];
                }
                cell_list[sum].Cell_fill(root); // заполняю ячейку
//-------------------------Нахождение шага сетки-----------------------------
                check_step(cell_list[sum].cell, this->max_diag);
//---------------------------------------------------------------------------
                sum++;                                                             // number of frames
            }
            beg_end[t][1] = sum - 1;
            t++;
        }
        beg_end_obj[i][1] = t - 1;
    }

    fin2 >> kl;                                                                     // число линий отрыва
    for (int i = 0; i < kl; i++) {
        fin2 >> kseg;                                                               // кол-во отрезков на каждой линии отрыва
        for (int j = 0; j < kseg; j++) {
            for (int k = 0; k < 2; k++) {
                fin2 >> x_bound[ssum][k][0] >> x_bound[ssum][k][1] >> x_bound[ssum][k][2];
            }
            ssum++;                                                                 //всего отрезков
        }
    }
    fin2.close();
    return;
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
    return;
}



