#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "element_geom.h"
#include "common_type_triangle.h"
#include "Cell_Triangle.h"



//===============================================================================================
//------------------------------------Constructor------------------------------------------------
//===============================================================================================
TGrid_DC_Full_Triangle::TGrid_DC_Full_Triangle(const std::string& filename) {
    
    grid_step = 0., num_frm = 0, num_mod = 0, num_blSeg = 0;

    int kmod, kseg, sign, nfm, not_used;
    int cur_frm = 0, cur_mod = 0;

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("Read geodat.dat error: Constructor");
    }

    // Заполнение сетки
    Cell_Triangle cell_tmp;
    double root[3][3];

    fin >> num_obj;
    for (int i = 0; i < num_obj; i++) {
        fin >> kmod;
        num_mod += kmod;
        std::array<int, 2> beo_i = {cur_mod, 0};

        for (int j = 0; j < kmod; j++) {
            fin >> sign >> nfm >> not_used >> not_used;
            num_frm += nfm;
            std::array<int, 2> be_i = {cur_frm, 0};
            for (int k = 0; k < nfm; k++) {
                for (int g = 0; g < 3; g++) {
                    fin >> root[g][0] >> root[g][1] >> root[g][2];
                }
                cell_tmp.cell_fill(root);
                cell_list.push_back(cell_tmp);
                check_step(cell_list[cur_frm].cell, this->grid_step);
                cur_frm++;
            }
            be_i[1] = cur_frm - 1;
            beg_endMod.push_back(be_i);
            cur_mod++;
        }
        beo_i[1] = cur_mod - 1;
        beg_endObj.push_back(beo_i);
    }


    // Заполнений линий отрыва
    fin >> num_bl;
    for (int i = 0; i < num_bl; i++) {
        fin >> kseg;
        num_blSeg += kseg;
        for (int j = 0; j < kseg; j++) {
            std::array<std::array<double, 3>, 2> xb_i;
            for (int k = 0; k < 2; k++) {
                fin >> xb_i[k][0] >> xb_i[k][1] >> xb_i[k][2];
            }
            x_bound.push_back(xb_i);
        }
    }
    fin.close();
}




//===============================================================================================
//-----------------------------------Copy Constructor--------------------------------------------
//===============================================================================================
TGrid_DC_Full_Triangle::TGrid_DC_Full_Triangle(const TGrid_DC_Full_Triangle& obj) {

    num_obj = obj.num_obj;
    num_mod = obj.num_mod;
    num_frm = obj.num_frm;
    grid_step = obj.grid_step;

    num_bl = obj.num_bl;
    num_blSeg = obj.num_blSeg;

    for (int i = 0; i < num_frm; i++) {
        cell_list.push_back(obj.cell_list[i]);
    }
   
    for (int i = 0; i < num_blSeg; i++) {
        x_bound.push_back(obj.x_bound[i]);
    }

    for (int i = 0; i < num_mod; i++) {
        beg_endMod.push_back(obj.beg_endMod[i]);
    }

    for (int i = 0; i < num_obj; i++) {
        beg_endObj.push_back(obj.beg_endObj[i]);
    }
}





//============================================================================================================
//------------------------------------------Destructor--------------------------------------------------------
//============================================================================================================
TGrid_DC_Full_Triangle::~TGrid_DC_Full_Triangle() {
    cell_list.clear();
    x_bound.clear();
    beg_endMod.clear();
    beg_endObj.clear();
}




//============================================================================================================
//------------------------------------------Step searching----------------------------------------------------
//============================================================================================================
void TGrid_DC_Full_Triangle::check_step(const double (&cur_root)[3][3], double cur_grid_step) {
    double diam_cur = get_diam(cur_root);
    if (diam_cur > cur_grid_step) {
        this->grid_step = diam_cur;
    }
}

